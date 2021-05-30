![SeismoVLAB Logo](../Logo.png)

Validations
===========

The validation folder provides with several models (from simple to sophisticated) that can be run using the **Pre-Analysis** and **Run-Analysis** modules. The validation folder is divided in the following sub-folders:

* `01-Debugging`: Folder that contains small simulation cases that can be run to test *Material*, *Element*, *Algorithm*, *Integrator* and *Solver* implementations. The **DEBUG CASE**'s names are as follows:
  **L01-Analysis_Formulation_Comment_Material_Element**, where: 
  * **L** is a letter that denotes complexity 
  * **Analysis** can be ST=Static, or DY=Dynamic
  * **Formulation** can be Lin=Linearized or Kin=kinematics
  * **Comment** is a short description 
  * **Material** and **Element** are the **SVL**'s class name.
* `02-Performance`: Folder that contains sophisticated simulation cases that can be run to test the parallel features of **SVL**.
* `03-Report`: Folder that run all debugging cases to verify the accuracy and well-behavior of the implemented features. 
  To generate the PDF report,  we advice to import the **01-Pre-Process** module as:

  ```bash
  export PYTHONPATH="/path/to/SVL/01-Pre_Process"
  ```

  Then execute in a python command prompt:

  ```python
  python3 '/path/to/runValidation.py'
  ```

All cases in folders `01-Debugging` and `02-Performance` are zipped (compressed); Therefore, they need to be unzipped before using them.

Further information can be obtained at:

* [Official website](https://www.seismovlab.com)
* [Download website](https://github.com/SeismoVLAB/SVL)
