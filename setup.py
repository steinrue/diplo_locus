import sys, os
from setuptools import setup, find_packages

conda_env_dir = os.environ.get('CONDA_PREFIX')
system_paths_dir = os.path.join(conda_env_dir, 'lib')
sys.path.append(system_paths_dir)

if __name__ == '__main__':
    setup(name='diplo_locus',
          description='Light-weight toolkit for the inference and simulation of Wright-Fisher diploid selection on independent loci from time-series data.',
          long_description=open('README.md', 'r').read(),
          long_description_content_type='text/markdown',
          version='1.1.0',
          python_requires=">=3.8",
          install_requires=['numpy', 'scipy', 'pandas', 'matplotlib', 'scikit-allel'],
          packages=find_packages(),  # 'src'
          py_modules=['diplo_locus', 'diplo_locus.likelihood', 'diplo_locus.utility',
                      'diplo_locus.simulate', 'diplo_locus.diffusion_core'],
          package_dir={'diplo_locus': 'src', 'src': 'src', '': 'src'},
          package_data={'': ['diplolocus_logo.png',
                             'src/README.md',
                             'examples/README.md',
                             'README.md']},
          install_package_data=True,
          scripts=['DiploLocus.py', 'DiploLocus_likelihood.py', 'DiploLocus_simulate.py'],
          entry_points={
              'console_scripts': ['DiploLocus = DiploLocus:main',
                                  'DiploLocus-likelihood = DiploLocus_likelihood:main',
                                  'DiploLocus-simulate = DiploLocus_simulate:main']
          },
          url='https://github.com/steinrue/diplo_locus',
          authors=['Xiaoheng Cheng', 'Matthias Steinruecken'],
          author_emails=['xhcheng@uchicago.edu', 'steinrue@uchicago.edu'],
          license='MIT'
          )
