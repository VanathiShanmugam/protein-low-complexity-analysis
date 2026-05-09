# Protein Low-Complexity Analysis

### Computational Extraction of Low-Complexity, Disorder, and Repeat Features from Protein Sequences

---

## Features

✔ Low-complexity region detection
✔ Disorder propensity estimation
✔ Repeat-density analysis
✔ FASTA-based feature extraction
✔ Batch processing of protein sequences
✔ CSV-based output generation

---

## Project Overview

This project presents a Python-based computational workflow for extracting sequence-derived protein features associated with:

* Low-complexity regions
* Intrinsic disorder propensity
* Repetitive sequence patterns

Protein sequences are parsed from FASTA files using Biopython, and multiple statistical features are computed for downstream bioinformatics and machine learning applications.

---

## ⚙️ Workflow

### 1. FASTA Parsing

Protein sequences are loaded using:

* Biopython `SeqIO`

---

### 2. Low-Complexity Analysis

Sliding-window analysis identifies regions dominated by repetitive amino acid composition.

Computed feature:

* `LowComplexity_Fraction`

---

### 3. Disorder Propensity Estimation

A FASTA-only disorder proxy estimates disorder-prone regions based on amino acid composition.

Computed features:

* `Mean_Disorder_Score`
* `Max_Disorder_Score`

---

### 4. Repeat Density Analysis

Repeated k-mer patterns are analyzed to estimate sequence repetitiveness.

Computed feature:

* `Repeat_Density`

---

### 5. Output Generation

Extracted features exported as:

```bash id="x2j9vq"
low_complexity_features.csv
```

---

## 🗂️ Project Structure

```bash id="l5d0ry"
protein-low-complexity-analysis/
│
├── README.md
├── requirements.txt
├── .gitignore
│
├── src/
│   └── extract_features.py
│
├── data/
│   └── sample.fasta
│
└── results/
```

---

## Usage

### Install dependencies

```bash id="m8q3af"
pip install -r requirements.txt
```

---

### Run analysis

```bash id="t9z5kh"
python src/extract_features.py
```

---

## 📁 Output

Generated CSV contains:

* Protein ID
* Low-complexity fraction
* Mean disorder score
* Maximum disorder score
* Repeat density

---

## Scientific Relevance

Low-complexity and disorder-associated regions are important in:

* Protein function prediction
* Pathogenicity analysis
* Intrinsically disordered protein studies
* Computational proteomics
* Feature engineering for ML pipelines

---

## Notes

* Uses simplified disorder estimation
* FASTA-only approach (no structural data required)
* Suitable for exploratory computational analysis

---

## Author

**Vanathi Shanmugam**
Bioinformatics | Genomics | Machine Learning

🔗 LinkedIn: https://www.linkedin.com/in/vanathi-shanmugam-26127928a

---

## 📜 License

Academic and research use only
