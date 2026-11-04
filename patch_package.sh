D=$(python -c "import snodiscov,os;print(os.path.dirname(snodiscov.__file__))")/workflow/envs
cat > $D/multiqc.yaml <<'EOF'
channels:
  - conda-forge
  - bioconda
  - nodefaults
dependencies:
  - python=3.11
  - setuptools<70
  - multiqc=1.14
EOF
sed -i 's/  - cutadapt=4.6/  - python=3.11\n  - cutadapt=4.6/' $D/cutadapt.yaml
