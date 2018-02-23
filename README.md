# HartreeFock
Very basic Hartree Fock implementation built to get a feeling for what happens under the hood. I wrote this in part as a
way to get my hands dirty with coding Hartree Fock from a computational/theoretical perspective, but also a way to learn 
Python. The code may therefore not be efficient or occassionally overloaded, just as a way for me to practice certain
Pythonic methods. 

I tried to write the code to be modular, so if you don't want to write code that parses basis set files to create a workable
molecule, then simply import the 'tools.py' file. The same goes for integrals, just import 'integrals.py'. Many people simply
prefer focusing on the scf part, in which case all you'll need to do is write your own 'scf.py' file.

To help the process, I've written up a document which gives everything you should need to implement your own Hartree Fock from
scratch. The file is aptly named, 'HartreeFockImplementationFromScratch.pdf'.

Prerequisites:
The code is written in Python 3, so please make sure you are using Python 3. That said, most of this can be easily changed to
fit Python 2.

Installing:
Simply 'git pull' the directory and everything should work, there is no installation necessary. All modules used should be 
readily available.
If this is not the case, you may need to 'pip install' certain modules.

General Structure:

'driver.py' -> Initializes the molecule using 'tools.py' and runs the scf procedure using 'scf.py'
'tools.py' -> Parses basis set files and coordinate files to create the molecule that you need in a program readable way
'scf.py' -> The meat of the program which builds all of the matrices (using 'integrals.py') and runs the actual scf routine
'integrals.py' -> Does the heavy lifting for the program, computing all of the intervals, from the simple overlap integrals
  to the more complicated eris.
'check.py' -> More or less a 'driver.py' file, but includes the ability to check against some test molecules

'test_molecules/' -> Houses the test molecules that give a range of test cases, atoms (Be), s-orbitals only (H2),
charged system (HeH+), p-orbital symmetric diatomic (N2), asymmetric diatomic (CO), tri-atomic (H2O), etc.
  
PLEASE let me know if I made mistakes anywhere. This whole project is a learning process and I'm sure I've made a number of 
mistakes along the way. I hope that it helps other interested students and future theoretical/compuational chemists and physicists get a better feel for what is means
to do ab-initio computation.
  
Best wishes for your learning!
Malte
langemf@uchicago.edu
