from setuptools import setup, find_packages
packages = find_packages(
    include=[
        'id_Dodd23', ], exclude=[
            '.txt', '.gitignore', '.npy', 'Old'])
setup(name='id_Dodd23',
      packages=packages
      )
