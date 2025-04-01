#!/bin/csh
source $ESP_ROOT/ESPenv.csh
ln -s $ESP_ROOT/externApps/Cart3Ddesign/demo/design/cart3d_demo.csm .
$ESP_ROOT/bin/ESPxddm model.input.xml /Model
if ($status != 0) then 
echo "ESPxddm Failed!"
exit 1
endif
exit 0
