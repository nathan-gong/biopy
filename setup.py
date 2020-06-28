import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="biopy-nathan-gong", 
    version="0.0.1",
    author="Nathan Gong",
    # author_email="nathangong9@gmail.com",
    description="A package containing common biology algorithms",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/nathan-gong/bio",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.X',
)