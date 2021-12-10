// ESP-python.js implements python functions for the Engineering Sketch Pad (ESP)
// written by John Dannenhoffer

// Copyright (C) 2010/2021  John F. Dannenhoffer, III (Syracuse University)
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

// interface functions that ESP-python.js provides
//    python.launch()
//    python.cmdUndo()
//    python.cmdSolve()
//    python.cmdSave()
//    python.cmdQuit()
//
//    python.cmdHome()                    not provided
//    python.cmdLeft()                    not provided
//    python.cmdRite()                    not provided
//    python.cmdBotm()                    not provided
//    python.cmdTop()                     not provided
//    python.cmdIn()                      not provided
//    python.cmdOut()                     not provided
//
//    python.mouseDown(e)                 not provided
//    python.mouseMove(e)                 not provided
//    python.mouseUp(e)                   not provided
//    python.mouseWheel(e)                not provided
//    python.mouseLeftCanvas(e)           not provided
//
//    python.keyPress(e)
//    python.keyDown(e)                   not provided
//    python.keyUp(e)                     not provided
//    python.keyPressPart1(myKeyPress)    not provided
//    python.keyPressPart2(picking, gprim)not provided
//    python.updateKeyWindow()
//
//    python.timLoadCB(text)              not provided
//    python.timSaveCB(text)              not provided
//    python.timQuitCB(text)
//    python.timMesgCB(text)

"use strict";


//
// callback when Python is launched
//
python.launch = function (e) {
    // alert("in python.launch(e="+e+")");

    // close the Tools menu
    var menu = document.getElementsByClassName("toolMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showToolMenu")) {
            menu[i].classList.remove(  "showToolMenu");
        }
    }

    // get the python file name
    var name;

    if (typeof e === "string") {
        name = e;
    } else {
        name = prompt("Enter name of python file", "../data/python/");
        if (name === null) {
            alert("A python file must be given");
            return
        }
    }

    // load the tim and run the python script
    browserToServer("timLoad|python|"+name+"|");

    // change mode
    changeMode(12);

    document.getElementById("toolMenuBtn").hidden = true;
    document.getElementById("doneMenuBtn").hidden = true;

    // change solve button legend
    var button = document.getElementById("solveButton");
    button["innerHTML"] = "RunPython";
    button.style.backgroundColor = "#3FFF3F";
};


//
// callback when "Python->Undo" is pressed (called by ESP.html)
//
python.cmdUndo = function () {
    alert("python.cmdUndo() is not implemented");
};


//
// callback when "solveButton" is pressed
//
python.cmdSolve = function () {
    // alert("in python.cmdSolve()");

    var button = document.getElementById("solveButton");
    button["innerHTML"] = "Running Python";
    button.style.backgroundColor = "#FFFF3F";

    browserToServer("timMesg|python|execute|")
};


//
// callback when "Python->Save" is pressed (called by ESP.html)
//
python.cmdSave = function () {
    alert("in python.cmdSave() is not implemented");
};


//
// callback when "Python->Quit" is pressed (called by ESP.html)
//
python.cmdQuit = function () {
    alert("in python.cmdQuit() is not implemented");
};


//
// callback when "homeButton" is pressed (calles by ESP.html)
//
//python.cmdHome = function () {
//    main.cmdHome();
//};


//
// callback when "leftButton" is pressed (calles by ESP.html)
//
//python.cmdLeft = function () {
//    main.cmdLeft();
//};


//
// callback when "riteButton" is pressed (calles by ESP.html)
//
//python.cmdRite = function () {
//    main.cmdRite();
//};


//
// callback when "botmButton" is pressed (calles by ESP.html)
//
//python.cmdBotm = function () {
//    main.cmdBotm();
//};


//
// callback when "topButton" is pressed (calles by ESP.html)
//
//python.cmdTop = function () {
//    main.cmdTop();
//};


//
// callback when "inButton" is pressed (calles by ESP.html)
//
//python.cmdIn = function () {
//    main.cmdIn();
//};


//
// callback when "outButton" is pressed (calles by ESP.html)
//
//python.cmdOut = function () {
//    main.cmdOut();
//};


//
// callback when any mouse is pressed in canvas
//
//python.mouseDown = function(e) {
//    main.mouseDown(e);
//};


//
// callback when any mouse moves in canvas
//
//python.mouseMove = function(e) {
//    main.mouseMove(e);
//};


//
// callback when the mouse is released in canvas
//
//python.mouseUp = function(e) {
//    main.mouseUp(e);
//};


//
// callback when the mouse wheel is rolled in canvas
//
//python.mouseWheel = function(e) {
//    main.mouseWheel(e);
//};


//
// callback when the mouse leaves the canvas
//
//python.mouseLeftCanvas = function(e) {
//    main.mouseLeftCanvas(e);
//};


//
// callback when a key is pressed
//
//python.keyPress = function (e) {
//    main.keyPress(e);
//};


//
// callback when an arror... or shift is pressed (needed for Chrome)
//
//python.keyDown = function (e) {
//    main.keyDown(e);
//};


//
// callback when a shift is released (needed for Chrome)
//
//python.keyUp = function (e) {
//    main.keyUp(e);
//};


//
// callback for first part of a keypress that is not recognized by wvUpdateUI
//
//python.keyPressPart1 = function(myKeyPress) {
//    // alert("in python.keyPressPart1(myKeyPress="+myKeyPress+")");
//    return 0;
//};


//
// callback for second part of a keypress that is not recognized by wvUpdateUI
//
//python.keyPressPart2 = function(picking, gprim) {
//    // alert("in python.keyPressPart2(picking="+picking+"   gprim="+gprim+")");
//};


//
// function to update the key window
//
//python.updateKeyWindow = function () {
//};


//
// callback when timLoad returns
//
python.timLoadCB = function (text) {
    // postMessage("in python.timLoadCB: "+text);
};


//
// callback when timSave returns
//
//python.timSaveCB = function (text) {
//    postMessage("in python.timSaveCB: "+text);
//};


//
// callback when timQuit returns
//
python.timQuitCB = function (text) {
    //postMessage("in plugs.timQuitCB: "+text);

    document.getElementById("toolMenuBtn").hidden = false;
    document.getElementById("doneMenuBtn").hidden = true;

    // change solve button legend
    var button = document.getElementById("solveButton");
    button["innerHTML"] = "Up to date";
    button.style.backgroundColor = "#FFFFFF";

    changeMode(0);
};


//
// callback when timMesg returns
//
python.timMesgCB = function (text) {
    // alert("in python.timMesgCB(text="+text+")")

    var tokens = text.split("|");

    // execute
    if        (tokens[1] == "execute") {

    // executeDone
    } else if (tokens[1] == "executeDone") {
        
        // update the display
        browserToServer("timDraw|python|");
        browserToServer("getBrchs|");
        browserToServer("getPmtrs|");

        // get stderr if it exists
        browserToServer("timMesg|python|stderr|");

        // set stdout if it exists
        browserToServer("timMesg|python|stdout|");

        // quit
        browserToServer("timQuit|python|");

        // change solve button legend
        var button = document.getElementById("solveButton");
        button["innerHTML"] = "Up to date";
        button.style.backgroundColor = "#FFFFFF";
        
    // stdout
    } else if (tokens[1] == "stdout") {
        if (tokens[2].length > 0) {
            postMessage( "Output of python program:\n"
                        +"*************************\n"
                        + text.substring(14)
                        +"*************************");
        } else {
            postMessage("Python program created no output");
        }

    // stderr
    } else if (tokens[1] == "stderr") {
        if (tokens[2].length > 0) {
            alert( "Error generated by python program:\n"
                   + text.substring(14, text.length-1));
        } else {
            postMessage("Python program did not raise an error");
        }

    // quit
    } else if (tokens[1] == "quit") {
        cmdQuit();

    // unknown message type
    } else {
        alert("unknown message type: "+ text);
    }
};

// /////////////////////////////////////////////////////////////////////

