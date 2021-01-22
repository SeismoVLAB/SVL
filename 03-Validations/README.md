![SeismoVLAB Logo](../Logo.png)

Validations
===========

The validation folder provides with several models (from simple to complicated) that can be run using the **Pre-Analysis** and **Run-Analysis** modules. The validation folder is divided in the following sub-folders:

* 01-Debugging: Several folder that contains small simulation cases that can be run to test *Material*, *Element*, *Algorithm*, *Integrator* and *Solver* implementations.  
* 02-Performance: Several folder that contains medium to large simulation cases that can be run to test the parallel features of **SVL**.
* 03-Report: In this folder, several debugging cases are presented in order to verify the accuracy and well-behaviour of the implemented features. The **DEBUG CASE**'s names are as follows:
  **L01-Analysis_Formulation_Comment_Material_Element**, where: 
  * **L** is a letter that denotes complexity 
  * **Analysis** can be ST=Static, or DY=Dynamic
  * **Formulation** can be Lin=Linearized or Kin=kinematics
  * **Comment** is a description 
  * **Material** and **Element** are the **SVL**'s class name.

All cases in folders **01-Debugging** and **02-Performance** are compressed; Therefore, they need to be unziped before using them.

To generate the PDF report,  we advice to import the 01-Pre-Process module as:

<pre>
export PYTHONPATH="/path/to/SeismoVLAB/01-Pre_Process"
</pre>

Then execute in a python command prompt:

<pre>
python3 runValidation.py
</pre>


Further information can be obatained at:

* Official website: https://www.seismovlab.com
* Download website: https://github.com/SeismoVLAB/SVL
