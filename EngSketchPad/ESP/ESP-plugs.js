// ESP-plugs.js implements plugs functions for the Engineering Sketch Pad (ESP)
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

// interface functions that ESP-plugs.js provides
//    plugs.launch()
//    plugs.cmdUndo()
//    plugs.cmdSolve()
//    plugs.cmdSave()
//    plugs.cmdQuit()
//
//    plugs.cmdHome()                    not provided
//    plugs.cmdLeft()                    not provided
//    plugs.cmdRite()                    not provided
//    plugs.cmdBotm()                    not provided
//    plugs.cmdTop()                     not provided
//    plugs.cmdIn()                      not provided
//    plugs.cmdOut()                     not provided
//
//    plugs.mouseDown(e)                 not provided
//    plugs.mouseMove(e)                 not provided
//    plugs.mouseUp(e)                   not provided
//    plugs.mouseWheel(e)                not provided
//    plugs.mouseLeftCanvas(e)           not provided
//
//    plugs.keyPress(e)                  not provided
//    plugs.keyDown(e)                   not provided
//    plugs.keyUp(e)                     not provided
//    plugs.keyPressPart1(myKeyPress)    not provided
//    plugs.keyPressPart2(picking, gprim)not provided
//    plugs.updateKeyWindow()
//
//    plugs.timLoadCB(text)
//    plugs.timSaveCB(text)
//    plugs.timQuitCB(text)
//    plugs.timMesgCB(text)

"use strict";


//
// name of TIM object
//
plugs.name = "plugs";


//
// callback when Plugs is launched
//
plugs.launch = function () {
    // alert("in plugs.launch()");

    // close the Tools menu
    var menu = document.getElementsByClassName("toolMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showToolMenu")) {
            menu[i].classList.remove(  "showToolMenu");
        }
    }

    // get name of cloudfile
    var name = prompt("Enter name of cloudfile");
    if (name === null) {
        alert("A cloudfile name must be given");
        return;
    }

    // load the tim
    browserToServer("timLoad|plugs|"+name+"|");

    // change done button legend
    var button = document.getElementById("doneMenuBtn");
    button["innerHTML"] = "Plugs";
    button.style.backgroundColor = "#3FFF3F";

    // change solve button legend
    button = document.getElementById("solveButton");
    button["innerHTML"] = "ExecutePhase1";
    button.style.backgroundColor = "#FFFF3F";

    // update the picture
    wv.sceneUpd = 1;

    // update the mode
    changeMode(11);
};


//
// callback when "Plugs->Undo" is pressed (called by ESP.html)
//
plugs.cmdUndo = function () {
    alert("plugs.cmdUndo() is not implemented");
};


//
// callback when "solveButton" is pressed
//
plugs.cmdSolve = function () {
    // alert("in plugs.cmdSolve()");

    // toggle between "ToggleToSize" and "ToggleToCent"
    var button = document.getElementById("solveButton");

    if (button["innerHTML"] == "ExecutePhase1") {
        plugs.npass = prompt("Enter number of Phase1 passes", "1");
        if (plugs.npass === null) {
            return;
        } else if (plugs.npass.search(/^[0-9]+$/i) < 0) {
            alert("npass (\""+plugs.npass+"\") must be an integer");
            return;
        } else if (plugs.npass == 0) {
            button["innerHTML"] = "ExecutePhase2"
            return;
        }

        plugs.ipass = 1;

        browserToServer("timMesg|plugs|phase1|");

        button["innerHTML"] =  "Executing pass "+plugs.ipass;
    } else if (button["innerHTML"] == "ExecutePhase2") {
        plugs.npass = prompt("Enter number of Phase2 passes", "25");
        if (plugs.npass === null) {
            return;
        } else if (plugs.npass.search(/^[0-9]+$/i) < 0) {
            alert("npass (\""+plugs.npass+"\") must be an integer");
            return;
        } else if (plugs.npass == 0) {
            button["innerHTML"] = "ExecutePhase1";
            return;
        }

        plugs.ipass = 1;

        browserToServer("timMesg|plugs|phase2|");

        button["innerHTML"] = "Executing pass "+plugs.ipass;
    } else {
        if (confirm("Plugs is currently executing.  Do you want to stop after this pass?") === true) {
            plugs.npass = 0;
        }
    }
};


//
// callback when "Plugs->Save" is pressed (called by ESP.html)
//
plugs.cmdSave = function () {
    // alert("in plugs.cmdSave()");

    // close the Plugs menu
    var menu = document.getElementsByClassName("doneMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showDoneMenu")) {
            menu[i].classList.remove(  "showDoneMenu");
        }
    }

    // execute save in the tim
    browserToServer("timSave|plugs|");

};


//
// callback when "Plugs->Quit" is pressed (called by ESP.html)
//
plugs.cmdQuit = function () {
    // alert("in plugs.cmdQuit()");

    // close the Plugs menu
    var menu = document.getElementsByClassName("doneMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showDoneMenu")) {
            menu[i].classList.remove(  "showDoneMenu");
        }
    }

    // execute quit in the tim
    browserToServer("timQuit|plugs|");

};


//
// callback when "homeButton" is pressed (calles by ESP.html)
//
//plugs.cmdHome = function () {
//    main.cmdHome();
//};


//
// callback when "leftButton" is pressed (calles by ESP.html)
//
//plugs.cmdLeft = function () {
//    main.cmdLeft();
//};


//
// callback when "riteButton" is pressed (calles by ESP.html)
//
//plugs.cmdRite = function () {
//    main.cmdRite();
//};


//
// callback when "botmButton" is pressed (calles by ESP.html)
//
//plugs.cmdBotm = function () {
//    main.cmdBotm();
//};


//
// callback when "topButton" is pressed (calles by ESP.html)
//
//plugs.cmdTop = function () {
//    main.cmdTop();
//};


//
// callback when "inButton" is pressed (calles by ESP.html)
//
//plugs.cmdIn = function () {
//    main.cmdIn();
//};


//
// callback when "outButton" is pressed (calles by ESP.html)
//
//plugs.cmdOut = function () {
//    main.cmdOut();
//};


//
// callback when any mouse is pressed in canvas
//
//plugs.mouseDown = function(e) {
//    main.mouseDown(e);
//}


//
// callback when any mouse moves in canvas
//
//plugs.mouseMove = function(e) {
//    main.mouseMove(e);
//}


//
// callback when the mouse is released in canvas
//
//plugs.mouseUp = function(e) {
//    main.mouseUp(e);
//}


//
// callback when the mouse wheel is rolled in canvas
//
//plugs.mouseWheel = function(e) {
//    main.mouseWheel(e);
//}


//
// callback when the mouse leaves the canvas
//
//plugs.mouseLeftCanvas = function(e) {
//    main.mouseLeftCanvas(e);
//}


//
// callback when a key is pressed
//
//plugs.keyPress = function (e) {
//    main.keyPress(e);
//}


//
// callback when an arror... or shift is pressed (needed for Chrome)
//
//plugs.keyDown = function (e) {
//    main.keyDown(e);
//};


//
// callback when a shift is released (needed for Chrome)
//
//plugs.keyUp = function (e) {
//    main.keyUp(e);
//};


//
// callback for first part of a keypress that is not recognized by wvUpdateUI
//
//plugs.keyPressPart1 = function(myKeyPress) {
//    // alert("in plugs.keyPressPart1(myKeyPress="+myKeyPress+")");
//    return 0;
//};


//
// callback for second part of a keypress that is not recognized by wvUpdateUI
//
//plugs.keyPressPart2 = function(picking, gprim) {
//    // alert("in plugs.keyPressPart2(picking="+picking+"   gprim="+gprim+")");
//};


//
// function to update the key window
//
plugs.updateKeyWindow = function () {
    var mesg = "Plugs\n\n";

    mesg += "classified   points " + plugs.classified;
    mesg += "unclassified points " + plugs.unclassified;
    mesg += "RMS error           " + plugs.RMS;

    var pre  = document.createElement("pre");
    var text = document.createTextNode(mesg);
    pre.appendChild(text);

    var plugsstat = document.getElementById("plugsStatus");

    plugsstat.replaceChild(pre, plugsstat.lastChild);
};


//
// callback when timLoad returns
//
plugs.timLoadCB = function (text) {
    // alert("in plugs.timLoadCB: "+text);

    var tokens = text.split('|');

    // if an error, get out of plugs
    if (tokens[2].substring(0,7) == "ERROR::") {
        alert("A problem occurred while trying to load \""+tokens[1]+"\"");
        plugs.cmdQuit();
    }
};


//
// callback when timSave returns
//
plugs.timSaveCB = function (text) {
    // alert("in plugs.timSaveCB: "+text);

    // change solve button legend
    var button = document.getElementById("solveButton");
    button["innerHTML"] = "Up to date";
    button.style.backgroundColor = null;

    changeMode(0);

    // rebuild
    browserToServer("build|-1|");
};


//
// callback when timQuit returns
//
plugs.timQuitCB = function (text) {
    // alert("in plugs.timQuitCB: "+text);

    // change solve button legend
    var button = document.getElementById("solveButton");
    button["innerHTML"] = "Up to date";
    button.style.backgroundColor = null;

    changeMode(0);

    // rebuild
    browserToServer("build|-1|");
};


//
// callback when timMesg returns
//
plugs.timMesgCB = function (text) {
    // alert("in plugs.timMesgCB: "+text);

    var button = document.getElementById("solveButton");
    var tokens = text.split("|");

    // phase1
    if (tokens[1] == "phase1") {

        if        (tokens[3] == 0) {
            postMessage("Phase1, pass "+plugs.ipass+" returned RMS="+tokens[2]);
            postMessage("Phase1 converged");
            plugs.npass = -1;
            button["innerHTML"] =  "ExecutePhase2";
        } else {
            postMessage("Phase1, pass "+plugs.ipass+" returned RMS="+tokens[2]+" (internal error "+tokens[3]+" detected)");
            return;
        }

        // update display
        browserToServer("timDraw|plugs|");

        browserToServer("timMesg|plugs|draw|");

        browserToServer("getPmtrs|");

        // if take another pass if there are more allowed
        plugs.ipass++;

        if (plugs.ipass <= plugs.npass) {
            sleep(100);

            button["innerHTML"] =  "Executing pass "+plugs.ipass;
            browserToServer("timMesg|plugs|phase1|");
        } else if (plugs.npass >= 0) {
            postMessage("Phase1 stopped because no passes available");
            button["innerHTML"] = "ExecutePhase1";
        }

    // phase2
    } else if (tokens[1] == "phase2") {

        var mesg = "Phase2, pass "+plugs.ipass+" returned unclass="+tokens[3]+", reclass="+tokens[4]+", RMS="+tokens[2];

        if        (tokens[5] == 0) {
            postMessage("Phase2 converged");
            plugs.npass = -1;
            button["innerHTML"] =  "ExecutePhase1";
        } else if (tokens[5] == 1) {
            postMessage(mesg+" (pass converged)");
        } else if (tokens[5] == 2) {
            postMessage(mesg+" (pass stalled)");
        } else if (tokens[5] == 3) {
            postMessage(mesg+" (pass ran out of iterations)");
        } else if (tokens[5] == 4) {
            postMessage(mesg+" (a DESPMTR has been perturbed)");
        } else if (tokens[5] == 5) {
            postMessage(mesg+" (no DESPMTR perturbation succeeded)");
            plugs.ipass = plugs.npass + 1;
            button["innerHTML"] =  "ExecutePhase1";
        } else {
            postMessage(mesg+" (internal error "+tokens[3]+" detected)");
            return;
        }

        // update display
        browserToServer("timDraw|plugs|");

        browserToServer("timMesg|plugs|draw|");

        browserToServer("getPmtrs|");

        // if take another pass if there are more allowed
        plugs.ipass++;

        if (plugs.ipass <= plugs.npass) {
            sleep(100);   // make sure display gets updated

            button["innerHTML"] =  "Executing pass "+plugs.ipass;
            browserToServer("timMesg|plugs|phase2|");
        } else if (plugs.npass >= 0) {
            postMessage("Phase2 stopped because no passes available");
            button["innerHTML"] = "ExecutePhase2";
        }

    // draw
    } else if (tokens[1] == "draw") {

    // unexpected message
    } else {
        postMessage("tokens[1]="+tokens[1]+" is an unexpected value");
        for (var i = 0; i < tokens.length; i++) {
            postMessage("tokens["+i+"]="+tokens[i]);
        }
    }
};

function sleep(miliseconds) {
   var currentTime = new Date().getTime();

   while (currentTime + miliseconds >= new Date().getTime()) {
   }
};


// /////////////////////////////////////////////////////////////////////

