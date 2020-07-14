from tixi3 import tixi3wrapper
from tigl3 import tigl3wrapper
import tigl3.configuration
import random
import sys

tixi = tixi3wrapper.Tixi3()
tigl = tigl3wrapper.Tigl3()
tixi.open('B7772VSP_v3.1.xml')
tigl.open(tixi,"")

mgr =  tigl3.configuration.CCPACSConfigurationManager_get_instance()

aircraft = mgr.get_configuration(tigl._handle.value)
fuselage = aircraft.get_fuselages().get_fuselage(1)
wings = aircraft.get_wings()

tigl.close()
tixi.close()
