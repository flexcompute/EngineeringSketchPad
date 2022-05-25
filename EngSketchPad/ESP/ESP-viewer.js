// ESP-viewer.js implements viewer functions for the Engineering Sketch Pad (ESP)
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

// interface functions that ESP-viewer.js provides
//    viewer.launch()
//    viewer.cmdUndo()
//    viewer.cmdSolve()
//    viewer.cmdSave()
//    viewer.cmdQuit()
//
//    viewer.cmdHome()                    not provided
//    viewer.cmdLeft()                    not provided
//    viewer.cmdRite()                    not provided
//    viewer.cmdBotm()                    not provided
//    viewer.cmdTop()                     not provided
//    viewer.cmdIn()                      not provided
//    viewer.cmdOut()                     not provided
//
//    viewer.mouseDown(e)                 not provided
//    viewer.mouseMove(e)                 not provided
//    viewer.mouseUp(e)                   not provided
//    viewer.mouseWheel(e)                not provided
//    viewer.mouseLeftCanvas(e)           not provided
//
//    viewer.keyPress(e)                  not provided
//    viewer.keyDown(e)                   not provided
//    viewer.keyUp(e)                     not provided
//    viewer.keyPressPart1(myKeyPress)    not provided
//    viewer.keyPressPart2(picking, gprim)not provided
//    viewer.updateKeyWindow()            not provided
//
//    viewer.timLoadCB(text)
//    viewer.timSaveCB(text)
//    viewer.timQuitCB(text)
//    viewer.timMesgCB(text)

"use strict";


//
// name of TIM object
//
viewer.name = "viewer";


//
// callback when Plugs is launched
//
viewer.launch = function () {
    alert("in viewer.launch()");

};


//
// callback when "Viewer->Undo" is pressed (called by ESP.html)
//
viewer.cmdUndo = function () {
    alert("viewer.cmdUndo() is not implemented");
};


//
// callback when "solveButton" is pressed
//
viewer.cmdSolve = function () {
    alert("viewer.cmdSolve() is not implemented");

};


//
// callback when "Viewer->Save" is pressed (called by ESP.html)
//
viewer.cmdSave = function () {
    alert("viewer.cmdSave() is not implemented");

};


//
// callback when "Viewer->Quit" is pressed (called by ESP.html)
//
viewer.cmdQuit = function () {
    alert("viewer.cmdQuit() is not implemented");

};


//
// callback when "homeButton" is pressed (calles by ESP.html)
//
viewer.cmdHome = function () {
    // alert("in viewer.cmdHome()");

    main.cmdHome();
};


//
// callback when "leftButton" is pressed (calles by ESP.html)
//
viewer.cmdLeft = function () {
    // alert("in viewer.cmdLeft()");

    main.cmdLeft();
};


//
// callback when "riteButton" is pressed (calles by ESP.html)
//
viewer.cmdRite = function () {
    // alert("in viewer.cmdRite()");

    main.cmdRite();
};


//
// callback when "botmButton" is pressed (calles by ESP.html)
//
viewer.cmdBotm = function () {
    // alert("in viewer.cmdBotm()");

    main.cmdBotm();
};


//
// callback when "topButton" is pressed (calles by ESP.html)
//
viewer.cmdTop = function () {
    // alert("in viewer.cmdTop()");

    main.cmdTop();
};


//
// callback when "inButton" is pressed (calles by ESP.html)
//
viewer.cmdIn = function () {
    // alert("in viewer.cmdIn()");

    main.cmdIn();
};


//
// callback when "outButton" is pressed (calles by ESP.html)
//
viewer.cmdOut = function () {
    // alert("in viewer.cmdOut()");

    main.cmdOut();
};


//
// callback when any mouse is pressed in canvas
//
viewer.mouseDown = function(e) {
    // alert("in viewer.mouseDown(e)");

    main.mouseDown(e);
}


//
// callback when any mouse moves in canvas
//
viewer.mouseMove = function(e) {
    // alert("in viewer.mouseMove(e)");

    main.mouseMove(e);
}


//
// callback when the mouse is released in canvas
//
viewer.mouseUp = function(e) {
    // alert("in viewer.mouseUp(e)");

    main.mouseUp(e);
}


//
// callback when the mouse wheel is rolled in canvas
//
viewer.mouseWheel = function(e) {
    // alert("in mouseWheel(e)");

    main.mouseWheel(e);
}


//
// callback when the mouse leaves the canvas
//
viewer.mouseLeftCanvas = function(e) {
    // alert("in mouseLeftCanvas(e)");

    main.mouseLeftCanvas(e);
}


//
// callback when a key is pressed
//
viewer.keyPress = function (keyPress) {
    // alert("in viewer.keyPress(keyPress="+keyPress+")");

    var myKeyPress = String.fromCharCode(keyPress.charCode);

    var done = 0;

    if        (myKeyPress == 'R' && wv.modifier == 1 && checkIfFree()) {
        browserToServer("timMesg|viewer|red|");
        done = 1;
    } else if (myKeyPress == 'G' && wv.modifier == 1 && checkIfFree()) {
        browserToServer("timMesg|viewer|green|");
        done = 1;
    } else if (myKeyPress == 'B' && wv.modifier == 1 && checkIfFree()) {
        browserToServer("timMesg|viewer|blue|");
        done = 1;
    }

    return done;
}


//
// callback when an arror... or shift is pressed (needed for Chrome)
//
viewer.keyDown = function (e) {
    // alert("in viewer.keyDown(e)");

    main.keyDown(e);
};


//
// callback when a shift is released (needed for Chrome)
//
viewer.keyUp = function (e) {
    // alert("in viewer.keyUp(e)");

    main.keyUp(e);
};


//
// callback for first part of a keypress that is not recognized by wvUpdateUI
//
viewer.keyPressPart1 = function(myKeyPress) {
    // alert("in viewer.keyPressPart1(myKeyPress="+myKeyPress+")");

    return 0;
};


//
// callback for second part of a keypress that is not recognized by wvUpdateUI
//
viewer.keyPressPart2 = function(picking, gprim) {
    // alert("in viewer.keyPressPart2(picking="+picking+"   gprim="+gprim+")");
};


//
// function to update the key window
//
//viewer.updateKeyWindow = function () {
//    // alert("in viewer.updateKeyWindow()");
//};


//
// callback when timQuit returns
//
viewer.timQuitCB = function (text) {
    // alert("in viewer.timQuitCB: "+text);

};


//
// callback when timMesg returns
//
viewer.timMesgCB = function (text) {
    // alert("in viewer.timMesgCB: "+text);

};

// /////////////////////////////////////////////////////////////////////

