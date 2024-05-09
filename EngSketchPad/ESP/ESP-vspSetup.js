// ESP-vspSetup.js implements vspSetup functions for the Engineering Sketch Pad (ESP)
// written by John Dannenhoffer

// Copyright (C) 2010/2024  John F. Dannenhoffer, III (Syracuse University)
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

// interface functions that ESP-vspSetup.js provides
//    vspSetup.launch()
//    vspSetup.cmdUndo()                    not provided
//    vspSetup.cmdSolve()                   not provided
//    vspSetup.cmdSave()                    not provided
//    vspSetup.cmdQuit()
//
//    vspSetup.cmdHome()                    not provided
//    vspSetup.cmdLeft()                    not provided
//    vspSetup.cmdRite()                    not provided
//    vspSetup.cmdBotm()                    not provided
//    vspSetup.cmdTop()                     not provided
//    vspSetup.cmdIn()                      not provided
//    vspSetup.cmdOut()                     not provided
//
//    vspSetup.mouseDown(e)                 not provided
//    vspSetup.mouseMove(e)                 not provided
//    vspSetup.mouseUp(e)                   not provided
//    vspSetup.mouseWheel(e)                not provided
//    vspSetup.mouseLeftCanvas(e)           not provided
//
//    vspSetup.keyPress(e)                  not provided
//    vspSetup.keyDown(e)                   not provided
//    vspSetup.keyUp(e)                     not provided
//    vspSetup.keyPressPart1(myKeyPress)    not provided
//    vspSetup.keyPressPart2(picking, gprim)not provided
//    vspSetup.updateKeyWindow()
//
//    vspSetup.timLoadCB(text)
//    vspSetup.timSaveCB(text)
//    vspSetup.timQuitCB(text)
//    vspSetup.timMesgCB(text)

"use strict";


//
// name of TIM object
//
vspSetup.name = "vspSetup";


//
// callback when VspSetup is launched
//
vspSetup.launch = function () {
    // alert("in vspSetup.launch()");

    // close the Tools menu
    var menu = document.getElementsByClassName("toolMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showToolMenu")) {
            menu[i].classList.remove(  "showToolMenu");
        }
    }

    // get name of vsp3file
    if (wv.filecomp.length > 0) {
        var defaultText = wv.filecomp;
        wv.filecomp = "";
    } else {
        defaultText = "";
    }
    
    vspSetup.vsp3name = prompt("Enter name of existing .vsp3 file (with * to complete)", defaultText);
    
    if (vspSetup.vsp3name === null) {
        alert("A vsp3 file name must be given");
        return;
    } else if (vspSetup.vsp3name.endsWith("*") === true) {
        browserToServer("completeFilename|vspSetup.launch()|"+vspSetup.vsp3name+"|");
        return;
    } else if (vspSetup.vsp3name.endsWith(".vsp3") === false) {
        vspSetup.vsp3name += ".vsp3";
    }

    // getUdcFile
    vspSetup.udcname = prompt("Enter name of the .udc file that you want to create/overwrite");
    if (vspSetup.udcname === null) {
        alert("A udc file name must be given");
        return;
    }
    if (vspSetup.udcname.endsWith(".udc") === false) {
        vspSetup.udcname += ".udc";
    }

    // load (and execute) the tim
    browserToServer("timLoad|vspSetup|"+vspSetup.vsp3name+"%"+vspSetup.udcname+"|");

    // since loading does all the work, quit immediately
    browserToServer("timQuit|vspSetup|");

    changeMode(16);
};


//
// callback when "VspSetup->Undo" is pressed (called by ESP.html)
//
vspSetup.cmdUndo = function () {
    alert("vspSetup.cmdUndo() is not implemented");
};


//
// callback when "solveButton" is pressed
//
vspSetup.cmdSolve = function () {
    alert("vspSetup.cmdSolve() is not implemented");
};


//
// callback when "VspSetup->Save" is pressed (called by ESP.html)
//
vspSetup.cmdSave = function () {
    // alert("in vspSetup.cmdSave()");

    // close the VspSetup menu
    var menu = document.getElementsByClassName("doneMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showDoneMenu")) {
            menu[i].classList.remove(  "showDoneMenu");
        }
    }

    // execute save in the tim
    browserToServer("timSave|vspSetup|");

};


//
// callback when "VspSetup->Quit" is pressed (called by ESP.html)
//
vspSetup.cmdQuit = function () {
    alert("vspSetup.cmdQuit() is not implemented");
};


//
// callback when "homeButton" is pressed (calles by ESP.html)
//
//vspSetup.cmdHome = function () {
//    main.cmdHome();
//};


//
// callback when "leftButton" is pressed (calles by ESP.html)
//
//vspSetup.cmdLeft = function () {
//    main.cmdLeft();
//};


//
// callback when "riteButton" is pressed (calles by ESP.html)
//
//vspSetup.cmdRite = function () {
//    main.cmdRite();
//};


//
// callback when "botmButton" is pressed (calles by ESP.html)
//
//vspSetup.cmdBotm = function () {
//    main.cmdBotm();
//};


//
// callback when "topButton" is pressed (calles by ESP.html)
//
//vspSetup.cmdTop = function () {
//    main.cmdTop();
//};


//
// callback when "inButton" is pressed (calles by ESP.html)
//
//vspSetup.cmdIn = function () {
//    main.cmdIn();
//};


//
// callback when "outButton" is pressed (calles by ESP.html)
//
//vspSetup.cmdOut = function () {
//    main.cmdOut();
//};


//
// callback when any mouse is pressed in canvas
//
//vspSetup.mouseDown = function(e) {
//    main.mouseDown(e);
//}


//
// callback when any mouse moves in canvas
//
//vspSetup.mouseMove = function(e) {
//    main.mouseMove(e);
//}


//
// callback when the mouse is released in canvas
//
//vspSetup.mouseUp = function(e) {
//    main.mouseUp(e);
//}


//
// callback when the mouse wheel is rolled in canvas
//
//vspSetup.mouseWheel = function(e) {
//    main.mouseWheel(e);
//}


//
// callback when the mouse leaves the canvas
//
//vspSetup.mouseLeftCanvas = function(e) {
//    main.mouseLeftCanvas(e);
//}


//
// callback when a key is pressed
//
//vspSetup.keyPress = function (e) {
//    main.keyPress(e);
//}


//
// callback when an arror... or shift is pressed (needed for Chrome)
//
//vspSetup.keyDown = function (e) {
//    main.keyDown(e);
//};


//
// callback when a shift is released (needed for Chrome)
//
//vspSetup.keyUp = function (e) {
//    main.keyUp(e);
//};


//
// callback for first part of a keypress that is not recognized by wvUpdateUI
//
//vspSetup.keyPressPart1 = function(myKeyPress) {
//    // alert("in vspSetup.keyPressPart1(myKeyPress="+myKeyPress+")");
//    return 0;
//};


//
// callback for second part of a keypress that is not recognized by wvUpdateUI
//
//vspSetup.keyPressPart2 = function(picking, gprim) {
//    // alert("in vspSetup.keyPressPart2(picking="+picking+"   gprim="+gprim+")");
//};


//
// function to update the key window
//
vspSetup.updateKeyWindow = function () {
    var mesg = "VspSetup\n\n";

    mesg += "classified   points " + vspSetup.classified;
    mesg += "unclassified points " + vspSetup.unclassified;
    mesg += "RMS error           " + vspSetup.RMS;

    var pre  = document.createElement("pre");
    var text = document.createTextNode(mesg);
    pre.appendChild(text);

    var vspSetupstat = document.getElementById("vspSetupStatus");

    vspSetupstat.replaceChild(pre, vspSetupstat.lastChild);
};


//
// callback when timLoad returns
//
vspSetup.timLoadCB = function (text) {
    // alert("in vspSetup.timLoadCB: "+text);

    var tokens = text.split('|');

    // if an error, get out of vspSetup
    if (tokens[2].substring(0,7) == "ERROR::") {
        alert("A problem occurred while trying to load \""+tokens[1]+"\"");
        vspSetup.cmdQuit();
    }
};


//
// callback when timSave returns
//
vspSetup.timSaveCB = function (text) {
    // alert("in vspSetup.timSaveCB: "+text);

    // change solve button legend
    var button = document.getElementById("solveButton");
    button["innerHTML"] = "Up to date";
    button.style.backgroundColor = null;

    changeMode(0);
};


//
// callback when timQuit returns
//
vspSetup.timQuitCB = function (text) {
    // alert("in vspSetup.timQuitCB: "+text);

    if (text.search(/ERROR::/) < 0) {

        // tell the user to add a UDPRIM call to the .csm file
        postMessage("To use the .udc file that was just created, add the following line to your .csm\n" +
                    "UDPRIM   /" + vspSetup.udcname.substring(0,vspSetup.udcname.length-4));
        alert("See message below to include \"" + vspSetup.vsp3name + "\" into your model");

    } else {
        alert("There was an error: "+text.substring(17));
    }

    // change solve button legend
    var button = document.getElementById("solveButton");
    button["innerHTML"] = "Up to date";
    button.style.backgroundColor = null;

    changeMode(0);
};


//
// callback when timMesg returns
//
vspSetup.timMesgCB = function (text) {
    alert("in vspSetup.timMesgCB: "+text);

//    var button = document.getElementById("solveButton");
    var tokens = text.split("|");

    // unexpected message
    postMessage("tokens[1]="+tokens[1]+" is an unexpected value");
    for (var i = 0; i < tokens.length; i++) {
        postMessage("tokens["+i+"]="+tokens[i]);
    }
};

// /////////////////////////////////////////////////////////////////////

