from setuptools import setup

if __name__ == "__main__":
    console_scripts = ["nmm = nmm:cli"]
    setup(
        cffi_modules="build_ext.py:ffibuilder",
        entry_points={"console_scripts": console_scripts},
    )
