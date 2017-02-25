Package: 

    * pbtranscript-internal-validation


Description: 

    * pbtranscript-internal-validation is a python package which contains 
      tools for validating E0, RC0 and human samples for PACBIO internal use ONLY.


INSTALLATION using pip:

    * `make pip-install`

    
INSTALLATION from pitchfork:

    * cd pitchfork
    * `make isoseq-core`
    * `make isoseq-validation`


INSTALLATION:

    * REQUIREMENT: pbtranscript must be installed
    * `python setup.py install --prefix=PATH_TO_YOUR_BUILD`


Usage:

    * `validate_smrtlink_isoseq_rc0.py path_to_a_smrtlink_isoseq_job validation_output_dir`
    * `validate_smrtlink_isoseq_rc0.py --help`


Document:

    * http://confluence.nanofluidics.com:8090/display/~yli/IsoSeq+Validation+and+Training

