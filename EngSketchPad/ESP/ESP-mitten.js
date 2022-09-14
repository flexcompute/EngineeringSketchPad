// ESP-mitten.js implements mitten functions for the Engineering Sketch Pad (ESP)
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

// interface functions that ESP-mitten.js provides
//    mitten.launch()
//    mitten.cmdUndo()
//    mitten.cmdSolve()
//    mitten.cmdSave()
//    mitten.cmdQuit()
//
//    mitten.cmdHome()                    not provided
//    mitten.cmdLeft()                    not provided
//    mitten.cmdRite()                    not provided
//    mitten.cmdBotm()                    not provided
//    mitten.cmdTop()                     not provided
//    mitten.cmdIn()                      not provided
//    mitten.cmdOut()                     not provided
//
//    mitten.mouseDown(e)                 not provided
//    mitten.mouseMove(e)                 not provided
//    mitten.mouseUp(e)                   not provided
//    mitten.mouseWheel(e)                not provided
//    mitten.mouseLeftCanvas(e)           not provided
//
//    mitten.keyPress(e)
//    mitten.keyDown(e)                   not provided
//    mitten.keyUp(e)                     not provided
//    mitten.keyPressPart1(myKeyPress)    not provided
//    mitten.keyPressPart2(picking, gprim)not provided
//    mitten.updateKeyWindow()
//
//    mitten.timLoadCB(text)              not provided
//    mitten.timSaveCB(text)              not provided
//    mitten.timQuitCB(text)              not provided
//    mitten.timMesgCB(text)

"use strict";


//
// name of TIM object
//
mitten.name = "mitten";


//
// callback when Mitten is launched
//
mitten.launch = function () {
    // alert("in mitten.launch()");

    // close the Tools menu
    var menu = document.getElementsByClassName("toolMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showToolMenu")) {
            menu[i].classList.remove(  "showToolMenu");
        }
    }

    // get component name
    var name = prompt("Enter name of Body");
    if (name === null) {
        alert("A Body name must be given");
        return;
    } else if (name.search(/^[a-zA-Z0-9]+$/i) < 0) {
        alert("Name (\""+name+"\") must only contain alphanumeric characters");
        return;
    }

    // load the tim
    browserToServer("timLoad|mitten|"+name+"|");
    
    // change done button legend
    var button = document.getElementById("doneMenuBtn");
    button["innerHTML"] = "Mitten";
    button.style.backgroundColor = "#3FFF3F";

    // change solve button legend
    button = document.getElementById("solveButton");
    button["innerHTML"] = "ToggleToSize";
    button.style.backgroundColor = "#FFFF3F";

    // update the picture
    wv.sceneUpd = 1;

    // update the mode
    changeMode(13);
};


//
// callback when "Mitten->Undo" is pressed (called by ESP.html)
//
mitten.cmdUndo = function () {
    alert("mitten.cmdUndo() is not implemented");
};


//
// callback when "solveButton" is pressed
//
mitten.cmdSolve = function () {
    // alert("in mitten.cmdSolve()");

    // toggle between "ToggleToSize" and "ToggleToCent"
    var button = document.getElementById("solveButton");
    if (button["innerHTML"] == "ToggleToSize") {
        button["innerHTML"] =  "ToggleToCent";

    } else {
        button["innerHTML"] = "ToggleToSize";
    }
    
    browserToServer(mesg+"0|");
};


//
// callback when "Mitten->Save" is pressed (called by ESP.html)
//
mitten.cmdSave = function () {
    // alert("in mitten.cmdSave()");

    // close the Mitten menu
    var menu = document.getElementsByClassName("doneMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showDoneMenu")) {
            menu[i].classList.remove(  "showDoneMenu");
        }
    }

    // execute save in the tim
    browserToServer("timSave|mitten|");

    // change solve button legend
    var button = document.getElementById("solveButton");
    button["innerHTML"] = "Up to date";
    button.style.backgroundColor = null;

    changeMode(0);

    // start the editor on the .csm file
    cmdFileEdit(null, 1);
};


//
// callback when "Mitten->Quit" is pressed (called by ESP.html)
//
mitten.cmdQuit = function () {
    // alert("in mitten.cmdQuit()");

    // close the Mitten menu
    var menu = document.getElementsByClassName("doneMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showDoneMenu")) {
            menu[i].classList.remove(  "showDoneMenu");
        }
    }

    // execute quit in the tim
    browserToServer("timQuit|mitten|");

    // change solve button legend
    var button = document.getElementById("solveButton");
    button["innerHTML"] = "Up to date";
    button.style.backgroundColor = null;

    changeMode(0);
};


//
// callback when "homeButton" is pressed (calles by ESP.html)
//
//mitten.cmdHome = function () {
//    main.cmdHome();
//};


//
// callback when "leftButton" is pressed (calles by ESP.html)
//
//mitten.cmdLeft = function () {
//    main.cmdLeft();
//};


//
// callback when "riteButton" is pressed (calles by ESP.html)
//
//mitten.cmdRite = function () {
//    main.cmdRite();
//};


//
// callback when "botmButton" is pressed (calles by ESP.html)
//
//mitten.cmdBotm = function () {
//    main.cmdBotm();
//};


//
// callback when "topButton" is pressed (calles by ESP.html)
//
//mitten.cmdTop = function () {
//    main.cmdTop();
//};


//
// callback when "inButton" is pressed (calles by ESP.html)
//
//mitten.cmdIn = function () {
//    main.cmdIn();
//};


//
// callback when "outButton" is pressed (calles by ESP.html)
//
//mitten.cmdOut = function () {
//    main.cmdOut();
//};


//
// callback when any mouse is pressed in canvas
//
//mitten.mouseDown = function(e) {
//    main.mouseDown(e);
//};


//
// callback when any mouse moves in canvas
//
//mitten.mouseMove = function(e) {
//    main.mouseMove(e);
//};


//
// callback when the mouse is released in canvas
//
//mitten.mouseUp = function(e) {
//    main.mouseUp(e);
//};


//
// callback when the mouse wheel is rolled in canvas
//
//mitten.mouseWheel = function(e) {
//    main.mouseWheel(e);
//};


//
// callback when the mouse leaves the canvas
//
//mitten.mouseLeftCanvas = function(e) {
//    main.mouseLeftCanvas(e);
//};


//
// callback when a key is pressed
//
mitten.keyPress = function (e) {
    // alert("in mitten.keyPress");

    var myKeyPress = String.fromCharCode(e.charCode);
    var button = document.getElementById("solveButton");

    postMessage("wv.keyPress="+wv.keyPress+"   myKeyPress="+myKeyPress);

    if        (myKeyPress == "x") {
        if (button["innerHTML"] == "ToggleToSize") {
            browserToServer("timMesg|mitten|xcent|-1|");
        } else {
            browserToServer("timMesg|mitten|xsize|0.5|");
        }

        return 1;
    } else if (myKeyPress == "X") {
        if (button["innerHTML"] == "ToggleToSize") {
            browserToServer("timMesg|mitten|xcent|+1|");
        } else {
            browserToServer("timMesg|mitten|xsize|2.0|");
        }

        return 1;
    } else if (myKeyPress == "y") {
        if (button["innerHTML"] == "ToggleToSize") {
            browserToServer("timMesg|mitten|ycent|-1|");
        } else {
            browserToServer("timMesg|mitten|ysize|0.5|");
        }

        return 1;
    } else if (myKeyPress == "Y") {
        if (button["innerHTML"] == "ToggleToSize") {
            browserToServer("timMesg|mitten|ycent|+1|");
        } else {
            browserToServer("timMesg|mitten|ysize|2.0|");
        }

        return 1;
    } else if (myKeyPress == "z") {
        if (button["innerHTML"] == "ToggleToSize") {
            browserToServer("timMesg|mitten|zcent|-1|");
        } else {
            browserToServer("timMesg|mitten|zsize|0.5|");
        }

        return 1;
    } else if (myKeyPress == "Z") {
        if (button["innerHTML"] == "ToggleToSize") {
            browserToServer("timMesg|mitten|zcent|+1|");
        } else {
            browserToServer("timMesg|mitten|zsize|2.0|");
        }

        return 1;
    } else if (myKeyPress == "R") {
        browserToServer("timMesg|mitten|rotate|0|");
        postMessage("rotating with angle = 0");

        return 1;
    } else if (myKeyPress == "C") {
        var seconds = prompt("Enter number of seconds", "10");
        if (seconds !== null) {
            browserToServer("timMesg|mitten|countdown|"+seconds+"|");
        }

        return 1;
    } else {
        postMessage("in mitten.keyPress(\""+myKeyPress+"\") passed back to main.keyPress");
    }

    // getting here means that we did not handle the keypress
    return 0;

};


//
// callback when an arror... or shift is pressed (needed for Chrome)
//
//mitten.keyDown = function (e) {
//    main.keyDown(e);
//};


//
// callback when a shift is released (needed for Chrome)
//
//mitten.keyUp = function (e) {
//    main.keyUp(e);
//};


//
// callback for first part of a keypress that is not recognized by wvUpdateUI
//
//mitten.keyPressPart1 = function(myKeyPress) {
//    // alert("in mitten.keyPressPart1(myKeyPress="+myKeyPress+")");
//    return 0;
//};


//
// callback for second part of a keypress that is not recognized by wvUpdateUI
//
//mitten.keyPressPart2 = function(picking, gprim) {
//    // alert("in mitten.keyPressPart2(picking="+picking+"   gprim="+gprim+")");
//};


//
// function to update the key window
//
mitten.updateKeyWindow = function () {
    var mesg = "Mitten\n\n";

    mesg +=    "xcent = " + mitten.xcent + "   xsize = " + mitten.xsize + "\n";

    mesg +=    "\nValid commands:\n";
    mesg +=    "  'x' 'y' 'z' to decrement\n";
    mesg +=    "  'X' 'Y' 'Z' to increment\n";

    var pre  = document.createElement("pre");
    var text = document.createTextNode(mesg);
    pre.appendChild(text);

    var mittenstat = document.getElementById("mittenStatus");
    
    mittenstat.replaceChild(pre, mittenstat.lastChild);
};


//
// callback when timLoad returns
//
//mitten.timLoadCB = function (text) {
//    postMessage("in mitten.timLoadCB: "+text);
//};


//
// callback when timSave returns
//
//mitten.timSaveCB = function (text) {
//    postMessage("in mitten.timSaveCB: "+text);
//};


//
// callback when timQuit returns
//
//mitten.timQuitCB = function (text) {
//    postMessage("in mitten.timQuitCB: "+text);
//};


//
// callback when timMesg returns
//
mitten.timMesgCB = function (text) {
    // postMessage("in mitten.timMesgCB: "+text);

    var tokens = text.split('|');
    if (tokens[1] == "rotate") {
        var angle  = Number(tokens[2]);
        if (angle < 360) {
            angle += 5;
            browserToServer("timMesg|mitten|rotate|"+angle+"|");
            postMessage("rotating with angle = "+angle);
        }
    }
};

// /////////////////////////////////////////////////////////////////////

