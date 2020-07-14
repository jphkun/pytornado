#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 10:46:00 2020

@author: cfse2
"""



from tixi3 import tixi3wrapper
from tigl3 import tigl3wrapper
from tigl3 import geometry

import tigl3.configuration



def print_point(p):
   print( "(" + str(p.x) + ";" + str(p.y) + ";" + str(p.z) + ")" )
  

def create_round_wing(wings, newWingUid, numberOfSection, diameter,sym):
  # parameters 
  deltaRotX = 180.0  / (numberOfSection - 1.0) 
  firstPosition = geometry.CTiglPoint(0,0,-diameter /2.0);
  firstNormal = geometry.CTiglPoint(0,-1,0);
  rot = geometry.CTiglTransformation();
  
  # create the wing
  wings.create_wing(newWingUid, numberOfSection , "NACA0012");
  wing = wings.get_wing(newWingUid);
  
  wing.set_symmetry(sym)

  
  # set the wing section elements  
  for idx  in range(1,wing.get_section_count() + 1) :
    rotX = (idx - 1) * deltaRotX;
    rot.set_identity(); 
    rot.add_rotation_x(rotX); 
    
    p = rot.transform(firstPosition);
    n = rot.transform(firstNormal); 
      
    s = wing.get_section(idx);
    e = s.get_section_element(1);
    ce = e.get_ctigl_section_element();
      
    ce.set_center(p);
    ce.set_normal(n); 
    if rotX >= 90:
      ce.set_rotation_around_normal(180); 

  
    
def save(tixi_h, aircraft, filename):
  aircraft.write_cpacs(aircraft.get_uid())
  configAsString = tixi_h.exportDocumentAsString();
  text_file = open(filename, "w")
  text_file.write(configAsString)
  text_file.close()
  
def create_space_ship():
    tixi_h = tixi3wrapper.Tixi3()
    tigl_h = tigl3wrapper.Tigl3()
    tixi_h.open("simpletest_cpacs_v3.1.xml")
    tigl_h.open(tixi_h, "")
 
    # get the configuration manager
    mgr =  tigl3.configuration.CCPACSConfigurationManager_get_instance()
    aircraft = mgr.get_configuration(tigl_h._handle.value)
    print(aircraft)
    wing0 = aircraft.get_wing(1);
    print(wing0)
    sym = wing0.get_symmetry();
    print(sym)
    wings = aircraft.get_wings(); 
    print(wings)
    
    
    create_round_wing(wings, "roundW", 29, 10,sym); 
    create_round_wing(wings, "roundW2", 13, 7,sym);
    create_round_wing(wings, "roundW3", 3, 8,sym);
    
    # wing3 = aircraft.get_wing("roundW3")
    # wing3.set_root_leposition(geometry.CTiglPoint(10,0,-4))
    
    # save(tixi_h, aircraft, "out-test.xml")
  
  
def main():
  create_space_ship(); 
  
if __name__ == '__main__':
  main(); 


# START: CFSE example
# from utils import *


# def create_round_wing(wings, newWingUid, numberOfSection, diameter,sym):
#     # parameters
#     deltaRotX = 180.0 / (numberOfSection - 1.0)
#     firstPosition = geometry.CTiglPoint(0,0,-diameter /2.0)
#     firstNormal = geometry.CTiglPoint(0,-1,0)
#     rot = geometry.CTiglTransformation()
  
#     # create the wing
#     wings.create_wing(newWingUid, numberOfSection, "NACA0012")
#     wing = wings.get_wing(newWingUid)
#     wing.set_symmetry(tigl3.core.TIGL_X_Z_PLANE)

#     # set the wing section elements  
#     for idx  in range(1,wing.get_section_count() + 1) :
#         rotX = (idx - 1) * deltaRotX;
#         rot.set_identity(); 
#         rot.add_rotation_x(rotX); 
    
#         p = rot.transform(firstPosition);
#         n = rot.transform(firstNormal); 
      
#         s = wing.get_section(idx);
#         e = s.get_section_element(1);
#         ce = e.get_ctigl_section_element();
      
#         ce.set_center(p);
#         ce.set_normal(n); 
#         if rotX >= 90:
#             ce.set_rotation_around_normal(180); 


# def main():
#     tigltixi = open_cpacs("B7772VSP_v3.1.xml")
#     tigl_h = tigltixi[0]
#     tixi_h = tigltixi[1]
#     aircraft = get_aircraft(tigl_h)
    
#     print_aircraft_info(aircraft)  
#     wing0 = aircraft.get_wing(1);
#     sym = wing0.get_symmetry();
#     # help(wing0)
#     # help(wing0.get_symmetry)
#     print("sym: " + str(sym) )
#     print("sym-type:" + str(type(sym)))
    
    
#     wings = aircraft.get_wings(); 
      
#     create_round_wing(wings, "roundW", 29, 10,sym); 
#     create_round_wing(wings, "roundW2", 13, 7,sym);
#     create_round_wing(wings, "roundW3", 3, 8,sym);
      
      
#     wing3 = aircraft.get_wing("roundW3")
#     wing3.set_root_leposition(geometry.CTiglPoint(10,0,-4))
      
#     save(tixi_h, aircraft, "out-test.xml")
      
      
#     tigl_h.close()
#     tixi_h.close()


# if __name__ == '__main__':
#     main()

# END: CFSE example

# # START: AIRFOIL EXMAPLE
# import tigl3.curve_factories
# import tigl3.surface_factories
# from OCC.gp import gp_Pnt
# from OCC.Display.SimpleGui import init_display
# from OCC.Geom import Handle_Geom_BSplineSurface
# import numpy as np

# # list of points on NACA2412 profile
# px = [1.000084, 0.975825, 0.905287, 0.795069, 0.655665, 0.500588, 0.34468, 0.203313, 0.091996, 0.022051, 0.0, 0.026892, 0.098987, 0.208902, 0.346303, 0.499412, 0.653352, 0.792716, 0.90373, 0.975232, 0.999916]
# py = [0.001257, 0.006231, 0.019752, 0.03826, 0.057302, 0.072381, 0.079198, 0.072947, 0.054325, 0.028152, 0.0, -0.023408, -0.037507, -0.042346, -0.039941, -0.033493, -0.0245, -0.015499, -0.008033, -0.003035, -0.001257]

# points_c1 = np.array([pnt for pnt in zip(px, [0.]*len(px), py)]) * 2.
# points_c2 = np.array([pnt for pnt in zip(px, [0]*len(px), py)])
# points_c3 = np.array([pnt for pnt in zip(px, py, [0.]*len(px))]) * 0.2

# # shift sections to their correct position
# # second curve at y = 7
# points_c2 += np.array([1.0, 7, 0])

# # third curve at y = 7.5
# points_c3[:, 1] *= -1
# points_c3 += np.array([1.7, 7.8, 1.0])

# curve1 = tigl3.curve_factories.interpolate_points(points_c1)
# curve2 = tigl3.curve_factories.interpolate_points(points_c2)
# curve3 = tigl3.curve_factories.interpolate_points(points_c3)

# surface = tigl3.surface_factories.interpolate_curves([curve1, curve2, curve3])
# # surface = tigl3.surface_factories.interpolate_curves([curve1, curve2, curve3], [0., 0.7, 1.])
# # surface = tigl3.surface_factories.interpolate_curves([curve1, curve2, curve3], degree=1)
# print(surface)
# # tigl3.surface_factories.interpolate_curves

# # start up the gui
# display, start_display, add_menu, add_function_to_menu = init_display()

# # make tesselation more accurate
# display.Context.SetDeviationCoefficient(0.0001)

# # draw the curve 
# # display.DisplayShape(curve1)
# # display.DisplayShape(curve2)
# # display.DisplayShape(curve3)
# display.DisplayShape(surface)

# # match content to screen and start the event loop
# display.FitAll()
# start_display()

# from ctypes import *
# from tixi3 import tixi3wrapper
# # define handles
# # tixiHandle = c_int(0)
# # xmlInputFilename    = "howtoin.xml"

# # open TIXI and TIGL shared libraries
# import sys
# TIXI = tixi3wrapper.Tixi3()

# Open a CPACS configuration file. First open the CPACS-XML file
# with TIXI to get a tixi handle.
# tixiReturn = TIXI.tixiOpenDocument(xmlInputFilename, byref(tixiHandle))
# if tixiReturn != 0:
#     print('Error: tixiOpenDocument failed for file: ' + xmlInputFilename)
#     exit(1)


# END: AIRFOIL EXMAPLE

# # -----------------------------------------------------------------------------
# # now demonstrate how to create a new document and adds some contents.
# # -----------------------------------------------------------------------------

# xmlOutputFilename = "howtoout.xml"
# handle = c_int(0)

# # Create a new document for writing with root element named
# # "plane". A file name is attributed to this document on saving.

# error = TIXI.create("plane", byref(handle))
# if error != 0:
#     # This error is fatal, exit here
#     quit()
# else:
#     print("Created a new xml file.")


# # First a header containing the tool name and version and the author
# # of the file is added. A timestamp is added automatically, as well.
# error = TIXI.tixiAddHeader(handle, "NoTool", "47.11", "Me")
# if error == 0:
#     print("Added header to new xml file.")

# # Now we insert a empty elements into the root element to create a
# # hierarchy of elements.
# nWings = 6
# error = TIXI.tixiAddTextElement(handle, "/plane", "wings", 0)
# if error == 0:
#     print("Inserted empty wings element.")

# # Insert empty wing elements into the wings element.
# for iWing in range(0, nWings):
#     error = TIXI.tixiAddTextElement(handle, "/plane/wings", "wing", 0)
#     if error == 0:
#         print("Inserted wing elements.")

# # Into each wing element insert a position attribute and a wingspan element.
# baseName = "/plane/wings/wing"

# # XPath indices start at 1 !!
# for iWing in range(1, nWings+1):
#     position = ""
#     if (iWing % 2) == 0:
#         position = "left"
#     else:
#         position = "right"

#     elementName = baseName + "[" + str(iWing) + "]"
#     # Add attribute
#     error = TIXI.tixiAddTextAttribute(handle,
#                                       elementName,
#                                       "position",
#                                       position)
#     if error == 0:
#         print("Added position attribute with value " + position + " to element " + elementName)

#     # Add wingspan element.
#     span = c_double(iWing + 0.1)
#     error = TIXI.tixiAddDoubleElement(handle,
#                                       elementName,
#                                       "wingspan",
#                                       span,
#                                       "%5.1f")
#     if error == 0:
#         print("Added wingspan element with value " + str(span.value) + " to element " + elementName)

# # write a 4x3 matrix
# nRows = 4
# nCols = 3
# dArray = c_double * (nRows * nCols)
# array = dArray()

# for i in range(0, nRows):
#     for j in range(0, nCols):
#         array[i*nCols+j] = i * 10. + j

# # Write array to matrix in row order (0)
# error = TIXI.tixiAddFloatMatrix(handle,
#                                 "/plane",
#                                 "four_by_three",
#                                 "rows",
#                                 "columns",
#                                 nRows,
#                                 nCols,
#                                 0,
#                                 array,
#                                 "%5.1f")
# if error == 0:
#     print("4x3 Matrix written")

# # Now write the array to the matrix in column order (1)
# error = TIXI.tixiAddFloatMatrix(handle,
#                                 "/plane",
#                                 "four_by_three",
#                                 "rows", "columns",
#                                 nRows,
#                                 nCols,
#                                 1,
#                                 array,
#                                 "%5.1f")
# if error == 0:
#     print("3x4 Matrix written")

# # Create an empty 3x3 matrix and fill the entries with a composite element.
# nRows = 3
# nCols = 3

# # first create an empty matrix
# error = TIXI.tixiCreateMatrix(handle,
#                               "/plane",
#                               "composite",
#                               "r",
#                               "c",
#                               nRows,
#                               nCols)

# if error == 0:
#     print("3x3 Matrix created.")

# # Path to a matrix entry is for example
# # /plane/composite/r[2]/c[3]
# basePath = "/plane/composite/r"

# for i in range(1, nRows + 1):
#     # build string for i-th the row element path
#     rowPath = basePath + "[" + str(i) + "]"

#     for j in range(1, nCols + 1):
#         position = ""
#         if (i*j) % 2 == 0:
#             position = "left"
#         else:
#             position = "right"

#         # build string for the path to matrix  entry (i,j)
#         entryName = rowPath + "/c[" + str(j) + "]"

#         # Add an empty wing element to the matrix element (i,j)
#         error = TIXI.tixiAddTextElement(handle, entryName, "wing", 0)
#         if error == 0:
#             print("entryName added.")

#         wingPath = entryName + "/wing"
#         # Insert a position and a span element into the wing element
#         error = TIXI.tixiAddTextElement(handle, wingPath, "position", position)
#         if error == 0:
#             print("wingPath/position added.")

#         span = c_double((i+j) * 10.)
#         error = TIXI.tixiAddDoubleElement(handle, wingPath, "span", span, 0)
#         if error == 0:
#             print("wingPath/span added.")

# # Now write a vector
# CArray = c_double * 10
# array = CArray(1,4,5.8,77.0,5,6,7,8,9,10)
# error = TIXI.tixiAddFloatVector(handle, "/plane", "testVector", array, 10)
# if error != 0:
#     print('Error in writing a vector with TIXI')
#     exit(1)

# # After all elements have been added, save the document to file
# error = TIXI.tixiSaveDocument(handle, xmlOutputFilename)
# if error == 0:
#     print("Document is written.")

# # Now we can close the document
# error = TIXI.tixiCloseDocument(handle)
# if error == 0:
#     print("TIXI document closed")
