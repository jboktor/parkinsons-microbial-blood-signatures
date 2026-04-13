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

---

## DeepRC._predict_proba wrong CSV separator (immuneML v3.0.21)

**File patched:** `immuneml_deeprc` conda env:
`$CONDA_PREFIX/lib/python3.11/site-packages/immuneML/ml_methods/classifiers/DeepRC.py`

**Line ~381 — `_predict_proba` method:**

Original:
```python
metadata_file_column_sep=DeepRCEncoder.SEP,
```

Patched:
```python
metadata_file_column_sep=DeepRCEncoder.METADATA_SEP,
```

**Why:** `DeepRCEncoder.SEP` is `"\t"` (the TSV separator for repertoire sequence files),
but the metadata file is a CSV (comma-separated). `_fit_for_label` correctly uses `","`,
but `_predict_proba` used `"\t"`, causing pandas to read the entire CSV row as a single
column. This made the `ID` column unfindable, producing:
```
KeyError: 'ID'
```

---

---

## DeepRC.make_data_loader None indices crash (immuneML v3.0.21)

**File patched:** `immuneml_deeprc` conda env:
`$CONDA_PREFIX/lib/python3.11/site-packages/immuneML/ml_methods/classifiers/DeepRC.py`

**Line ~214 — `make_data_loader` method:**

Added before `DeepRCRepDatasetSubset(...)`:
```python
if indices is None:
    indices = np.arange(len(full_dataset.target_features))
```

**Why:** `_predict_proba` calls `make_data_loader` with `indices=None` to use all samples,
but DeepRC's `DeepRCRepDatasetSubset.__init__` calls `np.asarray(indices, dtype=np.int)`
which fails on None. Must use `target_features` length (from metadata CSV), NOT
`len(full_dataset)` which returns `n_samples` from HDF5 (all samples, not just the subset).

---

---

## DeepRC._model_predict missing sequence_lengths arg (immuneML v3.0.21)

**File patched:** `immuneml_deeprc` conda env:
`$CONDA_PREFIX/lib/python3.11/site-packages/immuneML/ml_methods/classifiers/DeepRC.py`

**Line ~416 — `_model_predict` method:**

Original:
```python
logit_outputs = model(inputs, n_sequences)
```

Patched:
```python
logit_outputs = model(inputs, sequence_lengths, n_sequences)
```

**Why:** The `deeprc` package's `DeepRC.forward()` signature is
`forward(self, inputs_flat, sequence_lengths_flat, n_sequences_per_bag)` — 3 positional args.
The immuneML wrapper omitted `sequence_lengths`, causing:
```
TypeError: DeepRC.forward() missing 1 required positional argument: 'n_sequences_per_bag'
```
The training code in `deeprc.training.train` correctly passes all 3 args.

---

**Note:** All four patches must be reapplied if the `immuneml_deeprc` environment is recreated.
Consider reporting upstream to https://github.com/immuneML/immuneML/issues.
