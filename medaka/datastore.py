"""Storing of training and inference data to file."""
from collections import Counter, defaultdict, OrderedDict
from concurrent.futures import \
    as_completed, ProcessPoolExecutor, ThreadPoolExecutor
import warnings

import numpy as np
import yaml

import medaka.common

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    import h5py


class DataStore(object):
    """Read and write data to .hdf files."""

    _sample_path_ = 'samples'
    _groups_ = (
        'medaka_features_kwargs', 'medaka_model_kwargs', 'medaka_model_name',
        'medaka_label_decoding', 'medaka_feature_decoding',
        'medaka_label_counts', 'medaka_samples', 'medaka_multi_label',
        'medaka_label_description', 'medaka_label_scheme')

    def __init__(self, filename, mode='r', verify_on_close=True):
        """Initialize a datastore.

        :param filename: file to open.
        :param mode: file opening mode ('r', 'w', 'a').
        :param verify_on_close: on file close, check that all samples logged
            as being stored in file have a corresponding group within the
            `.hdf`."
        """
        self.filename = filename
        self.mode = mode
        self.verify_on_close = verify_on_close

        self._sample_keys = set()
        self.fh = None

        self.logger = medaka.common.get_named_logger('DataStore')

        self._meta = None
        self.write_executor = ThreadPoolExecutor(1)
        self.write_futures = []

    def __enter__(self):
        """Create filehandle."""
        self.fh = h5py.File(self.filename, self.mode)

        return self

    def __exit__(self, *args):
        """Verify file if requested."""
        if self.mode != 'r':
            if self.verify_on_close:
                self._verify_()  # verify data before saving meta
            else:
                self.logger.debug("Skipping validation on close.")
            self._write_metadata(self.meta)
            self.write_executor.shutdown(wait=True)
        self.fh.close()

    def _verify_(self):
        self.logger.debug("Verifying data.")
        self.fh.flush()
        fh = h5py.File(self.filename, 'r')
        # find the union of all present fields and remove and samples from the
        # index which don't contain all of these fields
        all_fields = set()

        for key in self.sample_keys:
            self.logger.debug("First round verify {}.".format(key))
            # if key is not in the file, remove it from the index
            grp = '{}/{}'.format(self._sample_path_, key)
            if grp not in fh:
                self.meta['medaka_samples'].remove(key)
                self.logger.debug(
                    "Removing sample {} as grp {} is not present.".format(
                        key, grp))
                continue
            all_fields.update(fh[grp].keys())

        for key in self.sample_keys:
            self.logger.debug("Second round verify {}.".format(key))
            for field in all_fields:
                path = '{}/{}/{}'.format(self._sample_path_, key, field)
                if path not in fh:
                    self.meta['medaka_samples'].remove(key)
                    self.logger.debug(
                        "Removing sample {} as field {} not present.".format(
                            key, path))
                    break

    @property
    def meta(self):
        """Meta data stored in file."""
        if self._meta is None:
            self._meta = self._load_metadata()
        return self._meta

    def update_meta(self, meta):
        """Update metadata."""
        self._meta = self.meta
        self._meta.update(meta)

    def write_sample(self, sample):
        """Write sample to hdf.

        Checks are performed to ensure a sample is not written twice and
        a count of unique training labels seen is maintained.

        :param sample: `medaka.common.Sample` object.
        """
        # count labels and store them in meta
        if 'medaka_label_counts' not in self.meta:
            self.meta['medaka_label_counts'] = Counter()
        # Store sample index in meta
        if 'medaka_samples' not in self.meta:
            self.meta['medaka_samples'] = set()

        if not any([
                isinstance(getattr(sample, field), np.ndarray)
                for field in sample._fields]):
            self.logger.debug('Not writing sample as it has no data.')
        elif sample.name not in self.meta['medaka_samples']:
            for field in sample._fields:
                if getattr(sample, field) is not None:
                    data = getattr(sample, field)
                    if isinstance(data, np.ndarray) and \
                            isinstance(data[0], np.unicode):
                        data = np.char.encode(data)
                    location = '{}/{}/{}'.format(
                        self._sample_path_, sample.name, field)
                    self.write_futures.append(
                        self.write_executor.submit(
                            self._write, location, data))

            if sample.labels is not None:
                # count combinations of labels accross haplotypes
                # cast nd.array to tuple of label tuples
                # use tolist to convert from np types to python types
                self.meta['medaka_label_counts'].update(
                    map(tuple, sample.labels.tolist()))
            # Do this last so we only add this sample to the index if we have
            # gotten this far
            self.meta['medaka_samples'].add(sample.name)
        else:
            self.logger.debug(
                'Not writing {} as present already'.format(
                    sample.name))

    def _write(self, location, data):
        # compress numpy arrays else write as is (scalars can't be compressed)
        if isinstance(data, np.ndarray):
            self.fh.create_dataset(
                location, data=data, compression='gzip', compression_opts=1)
        else:
            self.fh[location] = data

    def load_sample(self, key):
        """Load `medaka.common.Sample` object from file.

        :param key: str, sample name.
        :returns: `medaka.common.Sample` object.
        """
        s = {}
        for field in medaka.common.Sample._fields:
            pth = '{}/{}/{}'.format(self._sample_path_, key, field)
            if pth in self.fh:
                s[field] = self.fh[pth][()]
                if isinstance(s[field], np.ndarray) and \
                        isinstance(s[field][0], type(b'')):
                    s[field] = np.char.decode(s[field])
            else:
                s[field] = None
        return medaka.common.Sample(**s)

    def log_counts(self):
        """Log label counts."""
        def _label_name(label):
            if isinstance(label, tuple):
                return medaka.common.decoding[label[0]], label[1]
            else:
                return label
        self.logger.info("Label counts:\n{}".format('\n'.join(
            ['{}: {}'.format(_label_name(label), count)
                for label, count in self.meta['medaka_label_counts'].items()]
        )))

    def _write_metadata(self, data):
        """Save a data structure to file within a yml str."""
        self.logger.debug("Writing metadata.")
        for group, d in data.items():
            if group in self.fh:
                del self.fh[group]
            self.fh[group] = yaml.dump(d)

    def _load_metadata(self, groups=None):
        """Load meta data."""
        if groups is None:
            groups = self._groups_
        return {
            g: yaml.unsafe_load(self.fh[g][()])
            for g in groups if g in self.fh}

    @property
    def sample_keys(self):
        """Return tuple of sample keys."""
        value = tuple()
        if 'medaka_samples' in self.meta:
            value = tuple(self.meta['medaka_samples'])
        return value

    def _find_samples(self):
        """Find samples in a file and update meta with a list of samples."""
        samples = set()
        if self._sample_path_ in self.fh:
            samples = set(self.fh[self._sample_path_].keys())
        self.update_meta({'medaka_samples': samples})

    @property
    def n_samples(self):
        """Return the number of samples stored in file."""
        return len(self.sample_keys)


class DataIndex(object):
    """Index and serve samples from multiple `DataStore` compatible files."""

    def __init__(self, filenames, threads=4):
        """Intialize an index across a set of files.

        :param filenames: list of files to index.
        :param threads: number of threads to use for indexing.
        """
        self.logger = medaka.common.get_named_logger('DataIndex')

        self.filenames = filenames

        with DataStore(filenames[0]) as ds:
            self.logger.debug('Loading meta from {}'.format(filenames[0]))
            self.meta = ds.meta

        c_grp = 'medaka_label_counts'
        self.meta[c_grp] = Counter()

        if 'medaka_samples' in self.meta:
            del self.meta['medaka_samples']

        self.samples = []

        with ProcessPoolExecutor(threads) as executor:
            future_to_f = {executor.submit(
                DataIndex._load_meta, f): f for f in filenames}
            for i, future in enumerate(as_completed(future_to_f), 1):
                f = future_to_f[future]
                try:
                    meta = future.result()
                    if 'medaka_samples' in meta:
                        self.samples.extend(
                            [(s, f) for s in meta['medaka_samples']])
                        if c_grp in meta:
                            self.meta[c_grp].update(meta[c_grp])
                    else:
                        self.logger.info(
                            'Could not find samples in {}'.format(f))
                except Exception:
                    self.logger.info('Could not load meta from {}'.format(f))
                else:
                    self.logger.info(
                        'Loaded {}/{} ({:.2%}) sample files.'.format(
                            i, len(filenames), i / len(filenames)))

        # make order of samples independent of order in which tasks complete
        self.samples.sort()

        self._index = None

    @staticmethod
    def _load_meta(f):
        with DataStore(f) as ds:
            meta = ds.meta
            return meta

    @property
    def index(self):
        """Return a dictionary describing all samples.

        The dictionary maps sample names to their file and HDF group. It is
        sorted by reference coordinate.
        """
        if self._index is None:
            self._index = self._get_sorted_index()
        return self._index

    def _get_sorted_index(self):
        """Get index of samples indexed by reference and ordered by start pos.

        :returns: {ref_name: [sample dicts sorted by start]}

        """
        ref_names = defaultdict(list)

        for key, f in self.samples:
            d = medaka.common.Sample.decode_sample_name(key)
            if d is not None:
                d['key'] = key
                d['filename'] = f
                ref_names[d['ref_name']].append(d)

        # sort dicts so that refs are in order and within a ref,
        # chunks are in order
        ref_names_ordered = OrderedDict()

        # sort by start and -end so that if we have two samples with the same
        # start but differrent end points, the longest sample comes first
        def get_major_minor(x):
            return tuple((int(i) for i in x.split('.')))

        def sorter(x):
            return get_major_minor(
                x['start']) + tuple((-i for i in get_major_minor(x['end'])))

        for ref_name in sorted(ref_names.keys()):
            ref_names[ref_name].sort(key=sorter)
            ref_names_ordered[ref_name] = ref_names[ref_name]

        return ref_names_ordered

    def yield_from_feature_files(self, ref_names=None, samples=None):
        """Yield `medaka.common.Sample` objects from one or more feature files.

        :ref_names: iterable of str, only process these references.
        :samples: iterable of sample names to yield (in order in which they
            are supplied).

        :yields: `medaka.common.Sample` objects.

        """
        if samples is not None:
            # yield samples in the order they are asked for
            for sample, fname in samples:
                yield DataStore(fname).load_sample(sample)
        else:
            # yield samples sorted by ref_name and start
            if ref_names is None:
                ref_names = sorted(self.index.keys())
            for ref_name in ref_names:
                for d in self.index[ref_name]:
                    with DataStore(d['filename']) as ds:
                        yield ds.load_sample(d['key'])
