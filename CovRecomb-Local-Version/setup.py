'''
Author: error: git config user.name && git config user.email & please set dead value or install git
Date: 2022-07-04 14:33:43
LastEditors: Sonia-Ljy lijysunny@sina.com
LastEditTime: 2022-07-07 21:46:28
FilePath: /undefined/home/soniali/Desktop/03_recom_0531/0_code_org/CovRecomb-main/setup.py
Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
'''

from setuptools import setup, find_packages
from scripts import __version__, _program

setup(
    name='CovRecomb',
    version=__version__,
    packages=find_packages(),
    scripts=["scripts/function_set.py"],
    install_requires=["biopython>=1.70"],
    description='To identify the putative inter-lineage SARS-CoV-2 recombinants among consensus sequences inputted by users.',
    # url='https://github.com/wuaipinglab/CovRecomb-Local-Version.git',
    author='Jiaying Li',
    author_email='lijysunny@sina.com',
    entry_points="""
      [console_scripts]
      {program} = scripts.command:main
      """.format(program=_program),
    include_package_data=True,
    keywords=[],
    zip_safe=False)
