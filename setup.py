# coding=utf-8

from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()


setup(
    name = "cuteSVTrio",
    version = "0.1.0",
    description = "Long-read-based human genomic structural variation detection with cuteSVTrio",
    author = "Lixin",
    author_email = "xinli01@stu.hit.edu.cn",
    url = "https://github.com/tjiangHIT/cuteSVTrio",
    license = "MIT",
    packages = find_packages("src"),
    package_dir = {"": "src"},
    data_files = [("", ["LICENSE"])],
    scripts=['src/cuteSVTrio/cuteSVTrio'],
    # long_description = LONG_DESCRIPTION,
    long_description = readme,
    long_description_content_type = 'text/markdown',
    zip_safe = False,
    install_requires = ['scipy', 'pysam', 'Biopython', 'Cigar', 'numpy', 'pyvcf3', 'scikit-learn'],
)
