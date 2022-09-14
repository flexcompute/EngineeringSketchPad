// ESP-flowchart.js implements flowchart functions for the Engineering Sketch Pad (ESP)
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

// interface functions that ESP-flowchart.js provides
//    flowchart.launch()
//    flowchart.cmdUndo()
//    flowchart.cmdSolve()
//    flowchart.cmdSave()
//    flowchart.cmdQuit()
//
//    flowchart.cmdHome()                    not provided
//    flowchart.cmdLeft()                    not provided
//    flowchart.cmdRite()                    not provided
//    flowchart.cmdBotm()                    not provided
//    flowchart.cmdTop()                     not provided
//    flowchart.cmdIn()                      not provided
//    flowchart.cmdOut()                     not provided
//
//    flowchart.mouseDown(e)                 not provided
//    flowchart.mouseMove(e)                 not provided
//    flowchart.mouseUp(e)                   not provided
//    flowchart.mouseWheel(e)                not provided
//    flowchart.mouseLeftCanvas(e)           not provided
//
//    flowchart.keyPress(e)                  not provided
//    flowchart.keyDown(e)                   not provided
//    flowchart.keyUp(e)                     not provided
//    flowchart.keyPressPart1(myKeyPress)    not provided
//    flowchart.keyPressPart2(picking, gprim)not provided
//    flowchart.updateKeyWindow()
//
//    flowchart.timLoadCB(text)
//    flowchart.timSaveCB(text)
//    flowchart.timQuitCB(text)
//    flowchart.timMesgCB(text)

"use strict";


//
// name of TIM object
//
flowchart.name = "flowchart";


//
// callback when Flowchart is launched
//
flowchart.launch = function () {
    alert("in flowchart.launch()");

};


//
// callback when "Flowchart->Undo" is pressed (called by ESP.html)
//
flowchart.cmdUndo = function () {
    alert("flowchart.cmdUndo() is not implemented");
};


//
// callback when "solveButton" is pressed
//
flowchart.cmdSolve = function () {
    alert("in flowchart.cmdSolve() is not implemented");

};


//
// callback when "Flowchart->Save" is pressed (called by ESP.html)
//
flowchart.cmdSave = function () {
    alert("in flowchart.cmdSave() is not implemented");

};


//
// callback when "Flowchart->Quit" is pressed (called by ESP.html)
//
flowchart.cmdQuit = function () {
    alert("in flowchart.cmdQuit() is not implemented");

};


//
// callback when "homeButton" is pressed (calles by ESP.html)
//
//flowchart.cmdHome = function () {
//    main.cmdHome();
//};


//
// callback when "leftButton" is pressed (calles by ESP.html)
//
//flowchart.cmdLeft = function () {
//    main.cmdLeft();
//};


//
// callback when "riteButton" is pressed (calles by ESP.html)
//
//flowchart.cmdRite = function () {
//    main.cmdRite();
//};


//
// callback when "botmButton" is pressed (calles by ESP.html)
//
//flowchart.cmdBotm = function () {
//    main.cmdBotm();
//};


//
// callback when "topButton" is pressed (calles by ESP.html)
//
//flowchart.cmdTop = function () {
//    main.cmdTop();
//};


//
// callback when "inButton" is pressed (calles by ESP.html)
//
//flowchart.cmdIn = function () {
//    main.cmdIn();
//};


//
// callback when "outButton" is pressed (calles by ESP.html)
//
//flowchart.cmdOut = function () {
//    main.cmdOut();
//};


//
// callback when any mouse is pressed in canvas
//
//flowchart.mouseDown = function(e) {
//    main.mouseDown(e);
//}


//
// callback when any mouse moves in canvas
//
//flowchart.mouseMove = function(e) {
//    main.mouseMove(e);
//}


//
// callback when the mouse is released in canvas
//
//flowchart.mouseUp = function(e) {
//    main.mouseUp(e);
//}


//
// callback when the mouse wheel is rolled in canvas
//
//flowchart.mouseWheel = function(e) {
//    main.mouseWheel(e);
//}


//
// callback when the mouse leaves the canvas
//
//flowchart.mouseLeftCanvas = function(e) {
//    main.mouseLeftCanvas(e);
//}


//
// callback when a key is pressed
//
//flowchart.keyPress = function (e) {
//    main.keyPress(e);
//}


//
// callback when an arror... or shift is pressed (needed for Chrome)
//
//flowchart.keyDown = function (e) {
//    main.keyDown(e);
//};


//
// callback when a shift is released (needed for Chrome)
//
//flowchart.keyUp = function (e) {
//    main.keyUp(e);
//};


//
// callback for first part of a keypress that is not recognized by wvUpdateUI
//
//flowchart.keyPressPart1 = function(myKeyPress) {
//    // alert("in flowchart.keyPressPart1(myKeyPress="+myKeyPress+")");
//    return 0;
//};


//
// callback for second part of a keypress that is not recognized by wvUpdateUI
//
//flowchart.keyPressPart2 = function(picking, gprim) {
//    // alert("in flowchart.keyPressPart2(picking="+picking+"   gprim="+gprim+")");
//};


//
// function to update the key window
//
flowchart.updateKeyWindow = function () {
    alert("in flowchart.updateKeyWindow()");

};


//
// callback when timLoad returns
//
flowchart.timLoadCB = function (text) {
    alert("in flowchart.timLoadCB: "+text+" ignored)");

};


//
// callback when timSave returns
//
flowchart.timSaveCB = function (text) {
    alert("in flowchart.timSaveCB: "+text+" (ignored)");

};


//
// callback when timQuit returns
//
flowchart.timQuitCB = function (text) {
    alert("in flowchart.timQuitCB: "+text+" (ignored)");
};


//
// callback when timMesg returns
//
flowchart.timMesgCB = function (text) {
    //alert("in flowchart.timMesgCB: "+text);

    // open a new tab that contains the current flowchart
    window.open("../ESP/ESP-flowchart2.html");

};


// /////////////////////////////////////////////////////////////////////
