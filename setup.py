import setuptools


setuptools.setup(
    name="targetDB",
    version="0.0.1",
    author="Stephane De Cesco",
    author_email="stephane.decesco@ndm.ox.ac.uk",
    description="Package with an application to generate report on potential drug targets",
    long_description="Package with an application to generate report on potential drug targets (TO BE COMPLETED)",
    long_description_content_type="text/markdown",
    url="https://github.com/sdecesco/druggability",
    packages=setuptools.find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL3.0 License",
        "Operating System :: OS Independent",
    ),
	install_requires=[
		'biopython','scipy','matplotlib','pandas']
)
