// ESP-ereped.js implements ErepEd functions for the Engineering Sketch Pad (ESP)
// written by John Dannenhoffer

// Copyright (C) 2010/2022  John F. Dannenhoffer, III (Syracuse University)
//
// This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//    MA  02110-1301  USA

// interface functions that ESP-ereped.js provides
//    ereped.launch()
//    ereped.cmdUndo()
//    ereped.cmdSolve()
//    ereped.cmdSave()
//    ereped.cmdQuit()
//
//    ereped.cmdHome()                    not provided
//    ereped.cmdLeft()                    not provided
//    ereped.cmdRite()                    not provided
//    ereped.cmdBotm()                    not provided
//    ereped.cmdTop()                     not provided
//    ereped.cmdIn()                      not provided
//    ereped.cmdOut()                     not provided
//
//    ereped.mouseDown(e)                 not provided
//    ereped.mouseMove(e)                 not provided
//    ereped.mouseUp(e)                   not provided
//    ereped.mouseWheel(e)                not provided
//    ereped.mouseLeftCanvas(e)           not provided
//
//    ereped.keyPress(e)                  not provided
//    ereped.keyDown(e)                   not provided
//    ereped.keyUp(e)                     not provided
//    ereped.keyPressPart1(myKeyPress)
//    ereped.keyPressPart2(picking, gprim)
//
//    ereped.sceneUpdated()
//    ereped.updateKeyWindow()
//
//    ereped.timLoadCB(text)              not provided
//    ereped.timSaveCB(text)              not provided
//    ereped.timQuitCB(text)              not provided
//    ereped.timMesgCB(text)              not provided

// functions associated with ErepEd
//    erepInitialize(getNewName)
//    erepRedraw()

"use strict";


//
// name of TIM object
//
ereped.name = "ereped";


//
// callback when ErepEd is launched
//
ereped.launch = function () {
    // alert("in ereped.launch()");

    // close the Tools menu
    var menu = document.getElementsByClassName("toolMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showToolMenu")) {
            menu[i].classList.remove(  "showToolMenu");
        }
    }

    // make sure we do not enter if we are looking at an Erep already
    if (wv.plotType == 10) {
        alert("You cannot enter the Erep Editor while looking at an EBody\n"+
              "Switch the \"DisplayType\" back to \"0\" first");
        return;
    }

    // initialize GUI variables
    ereped.attrName  = "_erep";
    ereped.uncolored = 0;
    ereped.keptNodes = 0;
    ereped.keptEdges = 0;
    ereped.lastColor = -1;
    ereped.dihedral  = -1;
    ereped.colors    = [];

    // spectrum of colors
    ereped.spectrum = [0.87, 0.87, 0.87,       // light grey
                       1.00, 0.87, 0.87,       // light red
                       0.87, 1.00, 0.87,       // light green
                       0.87, 0.87, 1.00,       // light blue
                       0.87, 1.00, 1.00,       // light cyan
                       1.00, 0.87, 1.00,       // light magenta
                       1.00, 1.00, 0.87,       // light yellow
                       1.00, 0.50, 0.50,       // medium red
                       0.50, 1.00, 0.50,       // medium green
                       0.50, 0.50, 1.00,       // medium blue
                       0.50, 1.00, 1.00,       // medium cyan
                       1.00, 0.50, 1.00,       // medium magenta
                       1.00, 1.00, 0.50];      // medium yellow

    // uncolor all entities
    for (var gprim in wv.sceneGraph) {
        ereped.colors[gprim] = -1;
    }

    // initialize the application
    if (erepInitialize(1) > 0) {
        return;
    }

    // load the tim (which does nothing)
    browserToServer("timLoad|ereped|");

    // change done button legend
    var button = document.getElementById("doneMenuBtn");
    button["innerHTML"] = "ErepEd";
    button.style.backgroundColor = "#3FFF3F";

    // change solve button legend
    button = document.getElementById("solveButton");
    button["innerHTML"] = "ShowEBodys";
    button.style.backgroundColor = "#FFFF3F";

    // update the picture
    wv.sceneUpd = 1;

};


//
// callback when "ErepEd->Undo" is pressed (called by ESP.html)
//
ereped.cmdUndo = function () {
    alert("ereped.cmdUndo() is not implemented");
};


//
// callback when "solveButton" is pressed
//
ereped.cmdSolve = function () {
    // alert("in ereped.cmdSolve()");

    // if we are now showing the Ebody, switch back to showing the Bodys
    var button = document.getElementById("solveButton");
    if (button["innerHTML"] == "EditErep") {
        button["innerHTML"] = "ShowEBodys";

        // get the Body number
        var ibody = 0;
        for (var gprim in wv.sceneGraph) {
            var gprimList = gprim.trim().split(" ");

            if (gprimList[0] == "Body") {
                ibody = gprimList[1];
                break;
            }
        }

        // remove the Ebody
        browserToServer("timMesg|ereped|makeEBody|"+ibody+"|0|.|");

        // update the display
        browserToServer("timDraw|ereped|");

        return;
    }

    if (ereped.dihedral < 0) {
        var dihedral = prompt("Enter maximum dihedral angle", "5");
        if (dihedral === null) {
            return;
        } else if (isNaN(dihedral)) {
            alert("Illegal number format");
            return;
        } else if (dihedral < 0 || dihedral > 90) {
            alert("Out of range (0-90)");
            return;
        }

        ereped.dihedral = dihedral;
    }

    button["innerHTML"] = "EditErep";

    // get the Body number
    var ibody = 0;
    for (var gprim in wv.sceneGraph) {
        var gprimList = gprim.trim().split(" ");

        if (gprimList[0] == "Body") {
            ibody = gprimList[1];
            break;
        }
    }

    var mesg = "timMesg|ereped|makeEBody|" + ibody + "|" + ereped.dihedral + "|";


    // add Nodes to the mesg
    mesg += ereped.keptNodes + ";";

    for (var gprim in wv.sceneGraph) {
        var gprimList = gprim.trim().split(" ");

        if (gprimList[gprimList.length-2] == "Node") {
            if (ereped.colors[gprim] == -3) {
                mesg += gprimList[gprimList.length-1] + ";";
            }
        }
    }

    // add Edges to the mesg
    mesg += ereped.keptEdges + ";";

    for (var gprim in wv.sceneGraph) {
        var gprimList = gprim.trim().split(" ");

        if (gprimList[gprimList.length-2] == "Edge") {
            if (ereped.colors[gprim] == -2) {
                mesg += gprimList[gprimList.length-1] + ";";
            }
        }
    }

    // add Faces (in groups) to the mesg
    var colored = 0;
    for (gprim in wv.sceneGraph) {
        if (ereped.colors[gprim] > 0) {
            colored++;
        }
    }

    for (var igroup = 1; igroup < 99999; igroup++) {
        var nface = 0;

        for (gprim in wv.sceneGraph) {
            if (ereped.colors[gprim] == igroup) {
                nface++;
            }
        }

        if (nface > 1) {
            mesg += nface + ";";

            for (gprim in wv.sceneGraph) {
                if (ereped.colors[gprim] == igroup) {
                    gprimList = gprim.trim().split(" ");

                    mesg += gprimList[gprimList.length-1] + ";";
                    colored--;
                }
            }
        }

        if (colored <= 0) {
            break;
        }
    }

    // make the Ebody
    browserToServer(mesg+"0|");

    // update the display
    browserToServer("timDraw|ereped|");
};


//
// callback when "ErepEd->Save" is pressed (called by ESP.html)
//
ereped.cmdSave = function () {
    // alert("in ereped.cmdSave()");

    var ibrch    = brch.length;
    var nchange  = 0;
    var prevBody = -1;

    // close the ErepEd menu
    var menu = document.getElementsByClassName("doneMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showDoneMenu")) {
            menu[i].classList.remove(  "showDoneMenu");
        }
    }

    // we cannot save if we are currently viewing the EBodys
    var button = document.getElementById("solveButton");
    if (button["innerHTML"] != "ShowEBodys") {
        alert("You must be in editing mode to Save (press "+button["innerHTML"]+")");
        return;
    }

    // if we do not have a dihedral yet, ask for it now
    if (ereped.dihedral < 0) {
        var dihedral = prompt("Enter maximum dihedral angle", "5");
        if (dihedral === null) {
            return;
        } else if (isNaN(dihedral)) {
            alert("Illegal number format");
            return;
        } else if (dihedral < 0 || dihedral > 90) {
            alert("Out of range (0-90)");
            return;
        }

        ereped.dihedral = dihedral;
    }

    // send the changed attributes on Faces back to the server
    for (var gprim in wv.sceneGraph) {
        if (gprim.indexOf(" Face ") < 0) continue;

        var icolor = ereped.colors[gprim];
        if (icolor > 0) {

            // skip if we already have the attribute
            var attrs = wv.sgData[gprim];
            for (var i = 0; i < attrs.length; i+=2) {
                if (attrs[i] == ereped.attrName && attrs[i+1] == icolor) {
                    icolor = 0;
                    break;
                }
            }
            if (icolor == 0) {
                continue;
            }

            // add SELECT and ATTRIBUTE statements
            var nameList = gprim.trim().split(" ");

            if (nameList.length == 4) {
                postMessage("adding new atribute to "+gprim);

                if (nameList[1] != prevBody) {
                    browserToServer("newBrch|"+ibrch+"|select|$BODY|"+nameList[1]+"|||||||||");
                    ibrch++;

                    prevBody = nameList[1];

                    browserToServer("setAttr|"+ibrch+"|_erepAttr|1|$"+ereped.attrName+"|");
                    browserToServer("setAttr|"+ibrch+"|_erepAngle|1|"+ereped.dihedral+"|");
                }

                browserToServer("newBrch|"+ibrch+"|select|$FACE|"+nameList[3]+"|||||||||");
                ibrch++;

                browserToServer("setAttr|"+ibrch+"|"+ereped.attrName+"|1|"+icolor+"|");

                nchange++;
            } else {
                alert("could not set Attribute on named (unnumbered) Body");
            }
        }
    }

    // send the changes to the .Keep attribute on Edges back to the server
    for (var gprim in wv.sceneGraph) {
        if (gprim.indexOf(" Edge ") < 0) continue;

        try {
            var attrs = wv.sgData[gprim];

            // Edge is black (so it should not be kept)
            if (ereped.colors[gprim] == -1) {
                for (var i = 0; i < attrs.length; i+=2) {
                    if (attrs[i] == ".Keep") {
                        postMessage("removing .Keep from "+gprim);

                        var nameList = gprim.trim().split(" ");

                        browserToServer("newBrch|"+ibrch+"|select|$BODY|"+nameList[1]+"|||||||||");
                        ibrch++;

                        browserToServer("newBrch|"+ibrch+"|select|$EDGE|"+nameList[3]+"|||||||||");
                        ibrch++;

                        browserToServer("setAttr|"+ibrch+"|.Keep|1|<DeLeTe>|");

                        nchange++;

                        break;
                    }
                }

            // Edge is red (so it should be kept)
            } else {
                var found = 0;
                for (var i = 0; i < attrs.length; i+=2) {
                    if (attrs[i] == ".Keep") {
                        found++;
                        break;
                    }
                }

                if (found == 0) {
                    postMessage("adding .Keep to "+gprim);

                    var nameList = gprim.trim().split(" ");

                    browserToServer("newBrch|"+ibrch+"|select|$BODY|"+nameList[1]+"|||||||||");
                    ibrch++;

                    browserToServer("newBrch|"+ibrch+"|select|$EDGE|"+nameList[3]+"|||||||||");
                    ibrch++;

                    browserToServer("setAttr|"+ibrch+"|.Keep|1|1|");

                    nchange++;
                }
            }
        } catch (e) {
        }
    }

    // send the changes to the .Keep attribute on Nodes back to the server
    for (var gprim in wv.sceneGraph) {
        if (gprim.indexOf(" Node ") < 0) continue;

        // Node is black (so it should not be kept)
        if (ereped.colors[gprim] == -1) {
            for (var i = 0; i < attrs.length; i+=2) {
                if (attrs[i] == ".Keep") {
                    postMessage("removing .Keep from "+gprim);

                    var nameList = gprim.trim().split(" ");

                    browserToServer("newBrch|"+ibrch+"|select|$BODY|"+nameList[1]+"|||||||||");
                    ibrch++;

                    browserToServer("newBrch|"+ibrch+"|select|$NODE|"+nameList[3]+"|||||||||");
                    ibrch++;

                    browserToServer("setAttr|"+ibrch+"|.Keep|1|<DeLeTe>|");

                    nchange++;

                    break;
                }
            }

        // Node is red (so it should be kept)
        } else {
            var found = 0;
            for (var i = 0; i < attrs.length; i+=2) {
                if (attrs[i] == ".Keep") {
                    found++;
                    break;
                }
            }

            if (found == 0) {
                postMessage("adding .Keep to "+gprim);

                var nameList = gprim.trim().split(" ");

                browserToServer("newBrch|"+ibrch+"|select|$BODY|"+nameList[1]+"|||||||||");
                ibrch++;

                browserToServer("newBrch|"+ibrch+"|select|$NODE|"+nameList[3]+"|||||||||");
                ibrch++;

                browserToServer("setAttr|"+ibrch+"|.Keep|1|1|");

                nchange++;
            }
        }
    }

    // execute save in the tim (which does nothing)
    browserToServer("timSave|ereped|");

    // add a DUMP if any changes were made
//    if (nchange > 0) {
        browserToServer("newBrch|"+ibrch+"|dump|"+ereped.attrName+".egads|||||||||");

        postMessage("====> Re-build is needed <====");
//    }

    // update the key window if any changes were made
    if (nchange > 0) {
        ereped.updateKeyWindow();
    }

    changeMode(0);

    // get an updated version of the Branches and activate Build button
    wv.brchStat = 0;

    activateBuildButton();
};


//
// callback when "ErepEd->Quit" is pressed (called by ESP.html)
//
ereped.cmdQuit = function () {
    // alert("in ereped.cmdQuit()");

    // close the ErepEd menu
    var menu = document.getElementsByClassName("doneMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showDoneMenu")) {
            menu[i].classList.remove(  "showDoneMenu");
        }
    }

    // change solve button legend
    var button = document.getElementById("solveButton");
    button["innerHTML"] = "Up to date";
    button.style.backgroundColor = null;

    // execute quit in the tim (which does nothing)
    browserToServer("timQuit|ereped|");

    changeMode(0);

    // get original colors back
    browserToServer("build|-1|");
};


//
// callback when "homeButton" is pressed (calles by ESP.html)
//
//ereped.cmdHome = function () {
//    main.cmdHome();
//};


//
// callback when "leftButton" is pressed (calles by ESP.html)
//
//ereped.cmdLeft = function () {
//    main.cmdLeft();
//};


//
// callback when "riteButton" is pressed (calles by ESP.html)
//
//ereped.cmdRite = function () {
//    main.cmdRite();
//};


//
// callback when "botmButton" is pressed (calles by ESP.html)
//
//ereped.cmdBotm = function () {
//    main.cmdBotm();
//};


//
// callback when "topButton" is pressed (calles by ESP.html)
//
//ereped.cmdTop = function () {
//    main.cmdTop();
//};


//
// callback when "inButton" is pressed (calles by ESP.html)
//
//ereped.cmdIn = function () {
//    main.cmdIn();
//};


//
// callback when "outButton" is pressed (calles by ESP.html)
//
//ereped.cmdOut = function () {
//    main.cmdOut();
//};


//
// callback when any mouse is pressed in canvas
//
//ereped.mouseDown = function(e) {
//    main.mouseDown(e);
//}


//
// callback when any mouse moves in canvas
//
//ereped.mouseMove = function(e) {
//    main.mouseMove(e);
//}


//
// callback when the mouse is released in canvas
//
//ereped.mouseUp = function(e) {
//    main.mouseUp(e);
//}


//
// callback when the mouse wheel is rolled in canvas
//
//ereped.mouseWheel = function(e) {
//    main.mouseWheel(e);
//}


//
// callback when the mouse leaves the canvas
//
//ereped.mouseLeftCanvas = function(e) {
//    main.mouseLeftCanvas(e);
//}


//
// callback when a key is pressed
//
//ereped.keyPress = function (e) {
//    main.keyPress(e);
//};


//
// callback when an arror... or shift is pressed (needed for Chrome)
//
//ereped.keyDown = function (e) {
//    main.keyDown(e);
//};


//
// callback when a shift is released (needed for Chrome)
//
//ereped.keyUp = function (e) {
//    main.keyUp(e);
//};


//
// callback for first part of a keypress that is not recognized by wvUpdateUI
//
ereped.keyPressPart1 = function(myKeyPress) {
    // alert("in ereped.keyPressPart1(myKeyPress="+myKeyPress+")");
    console.trace();

    var button = document.getElementById("solveButton");
    var done = 0;

    if (myKeyPress == 'c') {
        if (button["innerHTML"] != "ShowEBodys") {
            alert("You must be in editing mode (press "+button["innerHTML"]+")");
            return;
        }

        wv.picking  = 99;
        wv.pick     = 1;
        wv.sceneUpd = 1;

        done = 1;
    } else if (myKeyPress == 'k') {
        if (button["innerHTML"] != "ShowEBodys") {
            alert("You must be in editing mode (press "+button["innerHTML"]+")");
            return;
        }

        wv.picking  = 107;
        wv.pick     = 1;
        wv.sceneUpd = 1;

        done = 1;
    } else if (myKeyPress == 'C') {
        if (button["innerHTML"] != "ShowEBodys") {
            alert("You must be in editing mode (press "+button["innerHTML"]+")");
            return;
        }

        var attrName = prompt("Enter Attribute name to import");
        if (attrName === null) {
            return;
        }

        var attrValu = prompt("Enter Attribute value to import (or -1 for all)");
        if (attrValu === null) {
            return;
        }

        // user provided an attribute value
        var nchange = 0;

        if (attrValu > 0) {
            var newColor = prompt("Enter <0 for first unused color\n"+
                                  "or    =0 to uncolor\n"+
                                  "or    >0 to select color", ereped.lastColor);
            if (newColor === null) {
                return;
            }

            // if user responded with negative, find first available color
            var icolor = 1;
            while (newColor < 0) {
                newColor = icolor;
                for (var gprim1 in wv.sceneGraph) {
                    if (ereped.colors[gprim1] == newColor) {
                        newColor = -1;
                        icolor++;
                        break;
                    }
                }
            }

            for (var gprim in wv.sceneGraph) {
                var gprimList = gprim.trim().split(" ");
                if        ((gprimList.length == 3 && gprimList[1] == "Face") ||
                           (gprimList.length == 4 && gprimList[2] == "Face")   ) {
                    try {
                        var attrs = wv.sgData[gprim];

                        for (var i = 0; i < attrs.length; i+=2) {
                            if (attrs[i] == attrName && Number(attrs[i+1]) == Number(attrValu)) {
                                if (ereped.colors[gprim] != newColor) {
                                    ereped.colors[gprim] = newColor;
                                    ereped.uncolored--;

                                    if (newColor == 0) {
                                        icolor = 0;
                                    } else {
                                        icolor = 1 + (newColor - 1) % 12;
                                    }
                                    wv.sceneGraph[gprim].fColor = [ereped.spectrum[3*icolor], ereped.spectrum[3*icolor+1], ereped.spectrum[3*icolor+2]];

                                    nchange++;
                                    postMessage("Changed color of "+gprim+" to "+newColor);
                                }
                            }
                        }
                    } catch (x) {
                        var attrs = [];
                    }
                }
            }

            // remember last color (for next time)
            ereped.lastColor = newColor;

        // user asked for all
        } else {
            var newColor = prompt("Enter first color (or -1 for fist available)", "-1");
            if (newColor === null) {
                return;
            }

            // if user responded with negative, find last color
            newColor = 0;
            for (var gprim1 in wv.sceneGraph) {
                if (ereped.colors[gprim1] > newColor) {
                    newColor = ereped.colors[gprim1];
                }
            }

            for (var gprim in wv.sceneGraph) {
                var gprimList = gprim.trim().split(" ");
                if        ((gprimList.length == 3 && gprimList[1] == "Face") ||
                           (gprimList.length == 4 && gprimList[2] == "Face")   ) {
                    try {
                        var attrs = wv.sgData[gprim];

                        for (var i = 0; i < attrs.length; i+=2) {
                            if (attrs[i] == attrName) {
                                var tempColor = Number(newColor) + Number(attrs[i+1]);
                                if (ereped.colors[gprim] != tempColor) {
                                    ereped.colors[gprim] = tempColor;
                                    ereped.uncolored--;

                                    icolor = 1 + (tempColor - 1) % 12;
                                    wv.sceneGraph[gprim].fColor = [ereped.spectrum[3*icolor], ereped.spectrum[3*icolor+1], ereped.spectrum[3*icolor+2]];

                                    nchange++;
                                    postMessage("Changed color of "+gprim+" to "+tempColor);
                                }
                            }
                        }
                    } catch (x) {
                        var attrs = [];
                    }
                }
            }
        }

        if (nchange > 0) {
            // update the Key Window
            ereped.updateKeyWindow();

            // update the picture
            wv.sceneUpd = 1;
        }

        done = 1;
    }

    return done;
};


//
// callback for second part of a keypress that is not recognized by wvUpdateUI
//
ereped.keyPressPart2 = function(picking, gprim) {
    // alert("in ereped.keyPressPart2(picking="+picking+"   gprim="+gprim+")");

    // second part of 'c' operation
    if (picking == 99) {
        if (gprim.indexOf(" Face ") < 0) {
            alert("Only Faces can be colored");
            return;
        }

        var oldColor  = ereped.colors[gprim];
        var nextColor = oldColor;

        // Face is currently not colored, so suggest last color (or 1 if first time)
        if (nextColor == 0) {
            if (ereped.lastColor < 0) {
                nextColor = 1;
            } else {
                nextColor = ereped.lastColor;
            }
        }

        // ask user for new color
        var newColor = prompt("Enter <0 for first unused color\n"+
                              "or    =0 to uncolor\n"+
                              "or    >0 to select color", nextColor);
        if (newColor === null) {
            return;
        }

        // if user responded with negative, find first available color
        var icolor = 1;
        while (newColor < 0) {
            newColor = icolor;
            for (var gprim1 in wv.sceneGraph) {
                if (ereped.colors[gprim1] == newColor) {
                    newColor = -1;
                    icolor++;
                    break;
                }
            }
        }

        // save the color and change the color of the Face
        ereped.colors[gprim] = newColor;
        postMessage("Changed color of "+gprim+" to "+newColor);

        if (newColor == 0) {
            icolor = 0;
        } else {
            icolor = 1 + (newColor - 1) % 12;
        }
        wv.sceneGraph[gprim].fColor = [ereped.spectrum[3*icolor], ereped.spectrum[3*icolor+1], ereped.spectrum[3*icolor+2]];

        // remember last color (for next time)
        ereped.lastColor = newColor;

        // adjust number of colored Faces
        if (oldColor == 0 && newColor > 0) {
            ereped.uncolored--;
        } else if (oldColor > 0 && newColor == 0) {
            ereped.uncolored++;
        }

        // update the Key Window
        ereped.updateKeyWindow();

    // second part of 'k' operation
    } else if (picking == 107) {
        if (gprim.indexOf(" Node ") < 0 && gprim.indexOf(" Edge ") < 0) {
            alert("Only Nodes and Edges can be kept");
            return;
        }

        var isEdge = 1;
        if (gprim.indexOf(" Node ") >= 0) {
            isEdge = 0;
        }

        // toggle the color of a Node
        if (isEdge == 0) {
            if (ereped.colors[gprim] == -1) {
                wv.sceneGraph[gprim].pColor = [1., 0., 0.];
                ereped.colors[gprim] = -3;

                postMessage("Marking "+gprim+" to be kept");
                ereped.keptNodes++;
                ereped.updateKeyWindow();
            } else {
                wv.sceneGraph[gprim].pColor = [0., 0., 0.];
                ereped.colors[gprim] = -1;

                postMessage("Unmarking "+gprim+" to be kept");
                ereped.keptNodes--;
                ereped.updateKeyWindow();
            }
        } else {
            if (ereped.colors[gprim] == -1) {
                wv.sceneGraph[gprim].lColor = [1., 0., 0.];
                ereped.colors[gprim] = -2;

                postMessage("Marking "+gprim+" to be kept");
                ereped.keptEdges++;
                ereped.updateKeyWindow();
            } else {
                wv.sceneGraph[gprim].lColor = [0., 0., 0.];
                ereped.colors[gprim] = -1;

                postMessage("Unmarking "+gprim+" to be kept");
                ereped.keptEdges--;
                ereped.updateKeyWindow();
            }
        }
    }
};


//
// callback when the scene graph has been updated
//
ereped.sceneUpdated = function () {
    // alert("in ereped.sceneUpdated()");

    var button = document.getElementById("solveButton");
    if (button["innerHTML"] != "EditErep") {
        erepRedraw();
    }
};


//
// function to update the key window
//
ereped.updateKeyWindow = function () {
    var mesg = "Erep Editor\n\n";

    mesg +=    "attrName          = " + ereped.attrName  + "\n";
    mesg +=    "# uncolored Faces = " + ereped.uncolored + "  (in grey)\n";
    mesg +=    "# kept Nodes      = " + ereped.keptNodes + "  (in red)\n";
    mesg +=    "# kept Edges      = " + ereped.keptEdges + "  (in red)\n\n";

    mesg +=    "Valid commands:\n";
    mesg +=    "  'c' color or uncolor a Face\n";
    mesg +=    "  'C' color Faces by attribute\n";
    mesg +=    "  'k' keep or unkeep a Node or Edge";

    var pre  = document.createElement("pre");
    var text = document.createTextNode(mesg);
    pre.appendChild(text);

    var timStatus = document.getElementById("timStatus");

    timStatus.replaceChild(pre, timStatus.lastChild);
}


//
// callback when timLoad returns
//
//ereped.timLoadCB = function (text) {
//    postMessage("in ereped.timLoadCB: "+text);
//}


//
// callback when timSave returns
//
//ereped.timSaveCB = function (text) {
//    postMessage("in ereped.timSaveCB: "+text);
//}


//
// callback when timQuit returns
//
//ereped.timQuitCB = function (text) {
//    postMessage("in ereped.timQuitCB: "+text);
//}


//
// callback when timMesg returns
//
//ereped.timMesgCB = function (text) {
//    postMessage("in ereped.timMesgCB: "+text);
//}

// /////////////////////////////////////////////////////////////////////


//
// Initialize ErepEd
//
var erepInitialize = function (getNewName) {
    // alert("in erepInitialize(getNewName="+getNewName+")");

    if (getNewName == 1) {
        var attrName = prompt("Enter attribute name:", ereped.attrName);
        if (attrName === null) {
            return 1;
        } else if (attrName.length <= 0) {
            return 1;
        } else {
            ereped.attrName = attrName;
        }
    }

    // make sure that any Faces that have ereped.attrName have it as
    //    only a single number
    for (var gprim in wv.sceneGraph) {
        var gprimList = gprim.trim().split(" ");
        if        ((gprimList.length == 3 && gprimList[1] == "Face") ||
                   (gprimList.length == 4 && gprimList[2] == "Face")   ) {
            if (ereped.colors[gprim] == -1) {
                ereped.colors[gprim] = 0;
                ereped.uncolored++;

                try {
                    var attrs = wv.sgData[gprim];

                    for (var i = 0; i < attrs.length; i+=2) {
                        if (attrs[i] == ereped.attrName) {
                            if (parseInt(attrs[i+1]) != Number(attrs[i+1])) {
                                alert("\""+ereped.attrName+"\" cannot be used because "+gprim+
                                      "\nalready has that attribute, and it is not a single integer");
                                ereped.attrname = "";
                                return 0;
                            } else if (Number(attrs[i+1]) < 0) {
                                alert("\""+ereped.attrName+"\" cannot be used because "+gprim+
                                      "\nalready has that attribute, and it has a negative value");
                                ereped.attrname = "";
                                return 0;
                            } else {
                                ereped.colors[gprim] = Number(attrs[i+1]);
                                ereped.uncolored--;
                            }
                        }
                    }
                } catch (e) {
                }
            }
        }
    }

    erepRedraw();

    // change mode
    changeMode(10);

    // update the Key Window
    ereped.updateKeyWindow();

    return 0;
};


//
// redraw the erep
//
var erepRedraw = function () {
    // alert("in erepRedraw()");

    for (var gprim in wv.sceneGraph) {
        var gprimList = gprim.trim().split(" ");

        // color Edges black (or red if they have a .Keep attribute */
        if ((gprimList.length == 3 && gprimList[1] == "Edge") ||
            (gprimList.length == 4 && gprimList[2] == "Edge")   ) {
            if (ereped.colors[gprim] == -2) {
                wv.sceneGraph[gprim].lColor = [1.0, 0.0, 0.0];
            } else {
                wv.sceneGraph[gprim].lColor = [0.0, 0.0, 0.0];
                ereped.colors[  gprim] = -1;

                try {
                    var attrs = wv.sgData[gprim];

                    for (var i = 0; i < attrs.length; i+=2) {
                        if (attrs[i] == ".Keep") {
                            ereped.keptEdges++;
                            wv.sceneGraph[gprim].lColor = [1.0, 0.0, 0.0];
                            ereped.colors[  gprim] = -2;
                            break;
                        }
                    }
                } catch (x) {
                }
            }

        // color Nodes black (or red if they have a .Keep attribute */
        } else if ((gprimList.length == 3 && gprimList[1] == "Node") ||
                   (gprimList.length == 4 && gprimList[2] == "Node")   ) {
            if (ereped.colors[gprim] == -3) {
                wv.sceneGraph[gprim].pColor = [1.0, 0.0, 0.0];
            } else {
                wv.sceneGraph[gprim].pColor = [0.0, 0.0, 0.0];
                ereped.colors[  gprim] = -1;

                try {
                    var attrs = wv.sgData[gprim];

                    for (var i = 0; i < attrs.length; i+=2) {
                        if (attrs[i] == ".Keep") {
                            ereped.keptNodes++;
                            wv.sceneGraph[gprim].pColor = [1.0, 0.0, 0.0];
                            ereped.colors[  gprim] = -3;
                            break;
                        }
                    }
                } catch (x) {
                }
            }
        }
    }

    // color the Faces
    for (gprim in wv.sceneGraph) {
        if (ereped.colors[gprim] == 0) {
            wv.sceneGraph[gprim].fColor = [ereped.spectrum[0], ereped.spectrum[1], ereped.spectrum[2]];
        } else if (ereped.colors[gprim] > 0) {
            var icolor = 1 + (ereped.colors[gprim]-1) % 12;
            wv.sceneGraph[gprim].fColor = [ereped.spectrum[3*icolor], ereped.spectrum[3*icolor+1], ereped.spectrum[3*icolor+2]];
        }
    }
};

// /////////////////////////////////////////////////////////////////////
