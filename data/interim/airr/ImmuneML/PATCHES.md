# ImmuneML Patches

## DeepRCEncoder ID column mismatch (immuneML v3.0.21)

**File patched:** `immuneml_deeprc` conda env:
`$CONDA_PREFIX/lib/python3.11/site-packages/immuneML/encodings/deeprc/DeepRCEncoder.py`

**Line 88 — `export_metadata_file` method:**

Original:
```python
metadata[DeepRCEncoder.ID_COLUMN] = dataset.get_repertoire_ids()
```

Patched:
```python
metadata[DeepRCEncoder.ID_COLUMN] = [f"{rid}.{DeepRCEncoder.EXTENSION}" for rid in dataset.get_repertoire_ids()]
```

**Why:** The HDF5 converter (`deeprc.dataset_converters.DatasetToHDF5`) stores sample keys
as `{uuid}.tsv` (matching the TSV filenames on disk), but `export_metadata_file` wrote bare
UUIDs without the `.tsv` extension into the metadata CSV `ID` column. When DeepRC's
`DeepRCRepDataset.__init__` tries to look up the metadata IDs in the HDF5 `sample_keys`,
they don't match, producing:
```
KeyError: "Samples ['' '' '' ...] could not be found in hdf5 file."
```

**Note:** This patch must be reapplied if the `immuneml_deeprc` environment is recreated.
Consider reporting upstream to https://github.com/immuneML/immuneML/issues.
