
from tixi3 import tixi3wrapper
from tigl3 import tigl3wrapper
from tigl3 import geometry

import tigl3.configuration

def open_cpacs(filename):
  
  tixi_h = tixi3wrapper.Tixi3()
  tixi_h.open(filename)
  value = tixi_h._handle.value
  
  tigl_h = tigl3wrapper.Tigl3()
 
  print("value" + str(value))
  print("version" + str(tixi_h.getDoubleElement("/cpacs/header/cpacsVersion")))
  tigl_h.open(tixi_h, "")
  return [tigl_h, tixi_h]
   
    
   
def save(tixi_h, aircraft, filename):
  aircraft.write_cpacs(aircraft.get_uid())
  configAsString = tixi_h.exportDocumentAsString();
  text_file = open(filename, "w")
  text_file.write(configAsString)
  text_file.close()
   
  
def get_aircraft(tigl_h):
  mgr =  tigl3.configuration.CCPACSConfigurationManager_get_instance()
  aircraft = mgr.get_configuration(tigl_h._handle.value)

  return aircraft  
  
  
def print_point(p):
   print( "(" + str(p.x) + ";" + str(p.y) + ";" + str(p.z) + ")" )
  
  
def print_aircraft_info(aircraft):
  print("aircraft: " + aircraft.get_uid())
  print("  wing count: " + str(aircraft.get_wing_count() ) );
  print("  fuselage count: " + str(aircraft.get_fuselage_count() ));
   
