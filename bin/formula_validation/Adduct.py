#!/usr/bin/env python

"""
This Python module contains not only the class Adduct, but also the test of
this Adduct class.

@contents :  This Python module contains not only the class Adduct, but also the test of
this Adduct class.
@project :  N/A
@program :  N/A
@file :  Adduct.py
@author :  Alberto Gil De la Fuente (alberto.gilf@gmail.com)

@version :  0.0.1, 20 July 2023
@information :  The Zen of Python
                  https://www.python.org/dev/peps/pep-0020/
                Style Guide for Python Code
                  https://www.python.org/dev/peps/pep-0008/
                Example NumPy Style Python Docstrings
                  http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html
                doctest – Testing through documentation
                  https://pymotw.com/2/doctest/

@copyright :  Copyright 2023 GNU AFFERO GENERAL PUBLIC.
              All rights are reserved. Reproduction in whole or in part is
              prohibited without the written consent of the copyright owner.
"""
__author__      = "Alberto Gil de la Fuente"
__copyright__   = "GPL License version 3"
from IncorrectAdduct import IncorrectAdduct
import re

class Adduct:

  """
  Methods:

  constructor(adduct): string representing an adduct in the form '[M+CH3CN+H]+', '[M-3H2O+2H]2+' or '[5M+Ca]2+' where the charge is specified at the end . It should start with a [, then contain the multimer number followed by an M, then the adduct formula with a +-, the closing ], and the number of charges indicated by a number and the symbol +-
  get_multimer(): number of multimers formed
  get_adduct_mass(): mass of the adduct
  get_adduct_charge(): number of charges, returns an integer > 0
  get_adduct_charge_type(): type of charge (+ or -)
  """ 
  __electron_weight=0.00054858

  def __init__(self, adduct: str):
    
    """
      Args:
        adduct (str): String representing an adduct in the form '[M+CH3CN+H]+', '[M-3H2O+2H]2+' or '[5M+Ca]2+' where the charge is specified at the end . It should start with a [, then contain the multimer number followed by an M, then the adduct formula with a +-, the closing ], and the number of charges indicated by a number and the symbol +-
      Returns:
        multimer (int): number of M in the structure where the adduct is present
        formula (str): formula or the adduct containing different subformulas starting with the symbol + or - to add or substract it
        charges: number of charges
        charge_type (+ or -): if the adduct is negatively or positively charged
      Raises:
        IncorrectAdduct: if the adduct does not follow the format specified
    """
    
    pattern = r'\[(\d)*M([\+-].*?)\](\d)*([\+-])'

    match = re.match(pattern, adduct)
    if match:
        self.__multimer = int(match.group(1).strip()) if match.group(1) is not None else 1
        self.__formula = match.group(2).strip()
        self.__adduct_mass = Adduct.__calculate_mass_from_formula_to_add_or_subtract(self.__formula)
        self.__charge = int(match.group(3).strip()) if match.group(3) is not None else 1
        self.__charge_type = match.group(4).strip() 
    else:
        raise IncorrectAdduct(adduct)
    
  def __eq__(self, other):
    if not isinstance(other, Adduct):
      return False
    return (self.__multimer == other.__multimer and self.__adduct_mass == other.__adduct_mass and self.__charge == other.__charge and self.__charge_type == other.__charge_type)

  def __str__(self):
    return f"Adduct(multimer={self.__multimer}, mass='{self.__adduct_mass}', charge={self.__charge_type}{self.__charge}, formula={self.__formula}"

  
  def __repr__(self):
    return str(self)

  def __hash__(self):
    return hash(self.__multimer, self.__adduct_mass, self.__charge, self.__charge_type)


  def __calculate_mass_from_formula_to_add_or_subtract(adduct_formula: str) -> float:
    """
      Args:
        formula (str): String representing the formula within an adduct in the form +HCOOH-H, +Ca, +H, +CH3COOH-H, etc..
      Returns:
        mass_from_adduct (int): the mass of the adduct to subtract to the structre
      Raises:
        IncorrectAdduct: if the adduct does not follow the format specified
    """
    from Formula import Formula
    mass_of_adduct_to_add_or_subtract=0
    pattern = r'([\+-])(\d)*([A-Za-z0-9]+)*'
    formulas = re.findall(pattern, adduct_formula)
    for symbol, number_subformulas, subformula_str in formulas:
      number_subformulas = int(number_subformulas) if number_subformulas else 1
      if symbol == '+':
        subformula=Formula.get_formula_from_str(subformula_str, 'None')
        monoisotopic_weight_subformula = subformula.get_monoisotopic_mass()
        mass_of_adduct_to_add_or_subtract += monoisotopic_weight_subformula * number_subformulas
      elif symbol == '-': 
        subformula=Formula.get_formula_from_str(subformula_str,'None')
        monoisotopic_weight_subformula = subformula.get_monoisotopic_mass()
        mass_of_adduct_to_add_or_subtract -= monoisotopic_weight_subformula * number_subformulas
      else: 
        raise IncorrectAdduct(adduct_formula)
      
    return mass_of_adduct_to_add_or_subtract
  
  def get_multimer(self) -> int:
    return self.__multimer

  def get_adduct_mass(self)-> float:
    return self.__adduct_mass
  
  def get_adduct_charge(self)-> int:
    return self.__charge

  def get_adduct_charge_type(self)-> str:
    return self.__charge_type

      
    

def main():
  import math
  print("=================================================================.")
  print("Test Case 1: Testing parsing from [M+H]+")
  print("=================================================================.")
  try:
    adduct_str1 = '[M+H]+'
    adduct_1 =Adduct(adduct_str1)
    if adduct_1.get_adduct_charge()==1 and adduct_1.get_adduct_charge_type() == '+' and adduct_1.get_multimer()==1 and math.isclose(adduct_1.get_adduct_mass(),1.0078,abs_tol=0.001):
      print("Test PASSED. The function to calculate the absolute values from ppms is correct.")
    else:
      print("Test FAILED. Check the function to calculate absolute values from ppms. Check adduct " + str(adduct_1) + " coming from: " + adduct_str1)
  except IncorrectAdduct as ia:
     print("Test FAILED. Check the function to calculate absolute values from ppms. Check adduct " + str(adduct_1) + " coming from: " + adduct_str1)
     print(str(ia))
  
  print("=================================================================.")
  print("Test Case 2: Testing parsing from [2M+HCOOH-H]-")
  print("=================================================================.")
  try:
    adduct_1 =Adduct('[2M+HCOOH-H]-')
    if adduct_1.get_adduct_charge()==1 and adduct_1.get_adduct_charge_type() == '-' and adduct_1.get_multimer()==2 and math.isclose(adduct_1.get_adduct_mass(),44.9976,abs_tol=0.001):
      print("Test PASSED. The function to calculate the absolute values from ppms is correct.")
    else:
      print("Test FAILED. Check the function to calculate absolute values from ppms. Check adduct " + str(adduct_1) + " coming from: " + adduct_str1)
  except IncorrectAdduct as ia:
     print("Test FAILED. Check the function to calculate absolute values from ppms. Check adduct " + str(adduct_1) + " coming from: " + adduct_str1)
     print(str(ia))

  print("=================================================================.")
  print("Test Case 3: [5M+Ca]2+")
  print("=================================================================.")
  try:
    adduct_1 =Adduct('[5M+Ca]2+')
    if adduct_1.get_adduct_charge()==2 and adduct_1.get_adduct_charge_type() == '+' and adduct_1.get_multimer()==5 and math.isclose(adduct_1.get_adduct_mass(),39.9625,abs_tol=0.001):
      print("Test PASSED. The function to calculate the absolute values from ppms is correct.")
    else:
      print("Test FAILED. Check the function to calculate absolute values from ppms. Check adduct " + str(adduct_1) + " coming from: " + adduct_str1)
  except IncorrectAdduct as ia:
     print("Test FAILED. Check the function to calculate absolute values from ppms. Check adduct " + str(adduct_1) + " coming from: " + adduct_str1)
     print(str(ia))

  print("=================================================================.")
  print("Test Case 4: [M-3H2O+2H]2+")
  print("=================================================================.")
  try:
    adduct_1 =Adduct('[M-3H2O+2H]2+')
    if adduct_1.get_adduct_charge()==2 and adduct_1.get_adduct_charge_type() == '+' and adduct_1.get_multimer()==1 and math.isclose(adduct_1.get_adduct_mass(),-52.016,abs_tol=0.001):
      print("Test PASSED. The function to calculate the absolute values from ppms is correct.")
    else:
      print("Test FAILED. Check the function to calculate absolute values from ppms. Check adduct " + str(adduct_1) + " coming from: " + adduct_str1)
  except IncorrectAdduct as ia:
     print("Test FAILED. Check the function to calculate absolute values from ppms. Check adduct " + str(adduct_1) + " coming from: " + adduct_str1)
     print(str(ia))


if __name__ == "__main__":
    main()