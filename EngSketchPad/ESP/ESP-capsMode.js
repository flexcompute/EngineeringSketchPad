// ESP-capsMode.js implements capsMode functions for the Engineering Sketch Pad (ESP)
// written by John Dannenhoffer

// Copyright (C) 2010/2025  John F. Dannenhoffer, III (Syracuse University)
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

// interface functions that ESP-capsMode.js provides
//    capsMode.launch()
//    capsMode.cmdUndo()
//    capsMode.cmdSolve()
//    capsMode.cmdSave()
//    capsMode.cmdQuit()
//
//    capsMode.cmdHome()                    not provided
//    capsMode.cmdLeft()                    not provided
//    capsMode.cmdRite()                    not provided
//    capsMode.cmdBotm()                    not provided
//    capsMode.cmdTop()                     not provided
//    capsMode.cmdIn()                      not provided
//    capsMode.cmdOut()                     not provided
//
//    capsMode.mouseDown(e)                 not provided
//    capsMode.mouseMove(e)                 not provided
//    capsMode.mouseUp(e)                   not provided
//    capsMode.mouseWheel(e)                not provided
//    capsMode.mouseLeftCanvas(e)           not provided
//
//    capsMode.keyPress(e)                  not provided
//    capsMode.keyDown(e)                   not provided
//    capsMode.keyUp(e)                     not provided
//    capsMode.keyPressPart1(myKeyPress)    not provided
//    capsMode.keyPressPart2(picking, gprim)not provided
//    capsMode.updateKeyWindow()
//
//    capsMode.timLoadCB(text)
//    capsMode.timSaveCB(text)
//    capsMode.timQuitCB(text)
//    capsMode.timMesgCB(text)
//
//    cmdCapsUpdateIntent()
//    cmdCapsCommitPhase()
//    cmdCapsQuitPhase()
//    cmdCapsSuspendPhase()
//    cmdCapsListPhases()
//    cmdCapsListHistory()
//    cmdCapsListAnalyses()
//    cmdCapsListBounds()
//    cmdCapsListAttachments()
//    cmdCapsAttachFile()
//    cmdCapsOpenAttachment()

"use strict";


//
// name of TIM object
//
capsMode.name = "capsMode";


//
// callback when CapsMode is launched
//
capsMode.launch = function () {
    // alert("in capsMode.launch()");

    // close the Tools menu
    var menu = document.getElementsByClassName("toolMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showToolMenu")) {
            menu[i].classList.remove(  "showToolMenu");
        }
    }

    // the following should never happer, but we want to make sure anyway
    if (wv.curTool != main || wv.curMode != 0) {
        alert("CapsMode cannot be entered when another mode ("+wv.curMode.name+") is active");
        return;
    }

    // get name of project:branch.revision
    if (wv.capsProj !== undefined && wv.capsPhase !== undefined) {
        var name = prompt("Enter project[:branch[.revision]][*]", wv.capsProj+":"+wv.capsPhase);
    } else {
        var name = prompt("Enter project[:branch[.revision]][*]");
    }
    if (name === null) {
        alert("A Project name must be given");
        return;
    } else if (name.length == 0) {
        alert("A Project name must be given");
        return;
    }

    // set the current tool so that capsMode.timLoadCB gets called
    //    if an error is encountered
    wv.curTool = capsMode;

    // load the tim
    capsMode.filenames = wv.filenames;
    var filenames = wv.filenames.split('|');
    browserToServer("timLoad|capsMode|"+name+"#"+filenames[1]+"|");

    // change done button legend and build the new menu
    var button = document.getElementById("doneMenuBtn");
    button["innerHTML"] = "Caps";
    button.style.backgroundColor = "#3FFF3F";

    var menu = document.getElementById("myDoneMenu");
    while (menu.firstChild != null) {
        menu.removeChild(menu.firstChild);
    }

    var button;

    button = document.createElement("input");
    button.type    = "button";
    button.title   = "Update Intent Phrase";
    button.value   = "Update Intent";
    button.onclick = cmdCapsUpdateIntent;
    menu.appendChild(button);

    menu.appendChild(document.createElement("hr"));

    button = document.createElement("input");
    button.type    = "button";
    button.title   = "List Phases";
    button.value   = "List Phases";
    button.onclick = cmdCapsListPhases;
    menu.appendChild(button);

    button = document.createElement("input");
    button.type    = "button";
    button.title   = "List Analyses";
    button.value   = "List Analyses";
    button.onclick = cmdCapsListAnalyses;
    menu.appendChild(button);

    button = document.createElement("input");
    button.type    = "button";
    button.title   = "List Bounds";
    button.value   = "List Bounds";
    button.onclick = cmdCapsListBounds;
    menu.appendChild(button);

    button = document.createElement("input");
    button.type    = "button";
    button.title   = "List History";
    button.value   = "List History";
    button.onclick = cmdCapsListHistory;
    menu.appendChild(button);

//    menu.appendChild(document.createElement("hr"));

//    button = document.createElement("input");
//    button.type    = "button";
//    button.title   = "Attach File";
//    button.value   = "Attach File";
//    button.onclick = cmdCapsAttachFile;
//    menu.appendChild(button);

//    button = document.createElement("input");
//    button.type    = "button";
//    button.title   = "List Attachments";
//    button.value   = "List Attachments";
//    button.onclick = cmdCapsListAttachments;
//    menu.appendChild(button);

//    button = document.createElement("input");
//    button.type    = "button";
//    button.title   = "Open Attachment";
//    button.value   = "Open Attachment";
//    button.onclick = cmdCapsOpenAttachment;
//    menu.appendChild(button);

    menu.appendChild(document.createElement("hr"));

    button = document.createElement("input");
    button.type    = "button";
    button.title   = "Commit Phase";
    button.value   = "Commit Phase";
    button.onclick = cmdCapsCommitPhase;
    menu.appendChild(button);

    button = document.createElement("input");
    button.type    = "button";
    button.title   = "Quit Phase";
    button.value   = "Quit Phase";
    button.onclick = cmdCapsQuitPhase;
    menu.appendChild(button);

    button = document.createElement("input");
    button.type    = "button";
    button.title   = "Suspend Phase";
    button.value   = "Suspend Phase";
    button.onclick = cmdCapsSuspendPhase;
    menu.appendChild(button);

    // change solve button legend and inactivate it
    button = document.getElementById("solveButton");
    button["innerHTML"] = "Loading Caps";
    button.style.backgroundColor = "#FFFF3F";

    // switch the server to CAPS mode and update the entries shown in the TreeWindow
    wv.server = "serveCAPS";
    rebuildTreeWindow(99);

    // update the mode
    changeMode(15);
};


//
// callback when "CapsMode->Undo" is pressed (called by ESP.html)
//
capsMode.cmdUndo = function () {
    alert("capsMode.cmdUndo() is not implemented");
};


//
// callback when "solveButton" is pressed
//
capsMode.cmdSolve = function () {
    // alert("capsMode.cmdSolve() is not implemented");
};


//
// callback when "CapsMode->Save" is pressed (called by ESP.html)
//
capsMode.cmdSave = function () {
    alert("in capsMode.cmdSave()");

};


//
// callback when "CapsMode->Quit" is pressed (called by ESP.html)
//
capsMode.cmdQuit = function () {
    alert("in capsMode.cmdQuit()");

};


//
// callback when "homeButton" is pressed (calles by ESP.html)
//
//capsMode.cmdHome = function () {
//    main.cmdHome();
//};


//
// callback when "leftButton" is pressed (calles by ESP.html)
//
//capsMode.cmdLeft = function () {
//    main.cmdLeft();
//};


//
// callback when "riteButton" is pressed (calles by ESP.html)
//
//capsMode.cmdRite = function () {
//    main.cmdRite();
//};


//
// callback when "botmButton" is pressed (calles by ESP.html)
//
//capsMode.cmdBotm = function () {
//    main.cmdBotm();
//};


//
// callback when "topButton" is pressed (calles by ESP.html)
//
//capsMode.cmdTop = function () {
//    main.cmdTop();
//};


//
// callback when "inButton" is pressed (calles by ESP.html)
//
//capsMode.cmdIn = function () {
//    main.cmdIn();
//};


//
// callback when "outButton" is pressed (calles by ESP.html)
//
//capsMode.cmdOut = function () {
//    main.cmdOut();
//};


//
// callback when any mouse is pressed in canvas
//
//capsMode.mouseDown = function(e) {
//    main.mouseDown(e);
//}


//
// callback when any mouse moves in canvas
//
//capsMode.mouseMove = function(e) {
//    main.mouseMove(e);
//}


//
// callback when the mouse is released in canvas
//
//capsMode.mouseUp = function(e) {
//    main.mouseUp(e);
//}


//
// callback when the mouse wheel is rolled in canvas
//
//capsMode.mouseWheel = function(e) {
//    main.mouseWheel(e);
//}


//
// callback when the mouse leaves the canvas
//
//capsMode.mouseLeftCanvas = function(e) {
//    main.mouseLeftCanvas(e);
//}


//
// callback when a key is pressed
//
//capsMode.keyPress = function (e) {
//    main.keyPress(e);
//}


//
// callback when an arror... or shift is pressed (needed for Chrome)
//
//capsMode.keyDown = function (e) {
//    main.keyDown(e);
//};


//
// callback when a shift is released (needed for Chrome)
//
//capsMode.keyUp = function (e) {
//    main.keyUp(e);
//};


//
// callback for first part of a keypress that is not recognized by wvUpdateUI
//
//capsMode.keyPressPart1 = function(myKeyPress) {
//    // alert("in capsMode.keyPressPart1(myKeyPress="+myKeyPress+")");
//    return 0;
//};


//
// callback for second part of a keypress that is not recognized by wvUpdateUI
//
//capsMode.keyPressPart2 = function(picking, gprim) {
//    // alert("in capsMode.keyPressPart2(picking="+picking+"   gprim="+gprim+")");
//};


//
// function to update the key window
//
capsMode.updateKeyWindow = function () {
    var mesg = "CapsMode\n\n";

    mesg += "classified   points " + capsMode.classified;
    mesg += "unclassified points " + capsMode.unclassified;
    mesg += "RMS error           " + capsMode.RMS;

    var pre  = document.createElement("pre");
    var text = document.createTextNode(mesg);
    pre.appendChild(text);

    var capsModestat = document.getElementById("capsModeStatus");

    capsModestat.replaceChild(pre, capsModestat.lastChild);
};


//
// callback when timLoad returns
//
capsMode.timLoadCB = function (text) {
    // alert("in capsMode.timLoadCB: "+text);

    var tokens = text.split('|');

    // check to see if the requested Phase is locked
    if (tokens[1].substring(0,7) == "ERROR::" && tokens[1].includes(" is locked by:")) {
        if (confirm(tokens[1].substring(8)) !== true) {
            alert("Cannot continue because the Phase is locked");

            // return the done menu to its original state
            var button = document.getElementById("doneMenuBtn");
            button["innerHTML"] = "Done";
            button.style.backgroundColor = null;

            var menu = document.getElementById("myDoneMenu");
            while (menu.firstChild != null) {
                menu.removeChild(menu.firstChild);
            }

            // return the solve button to its original state
            button = document.getElementById("solveButton");
            button["innerHTML"] = "Up to date";
            button.style.backgroundColor = null;

            // return to "main" mode (which is the only mode from which
            //    this TIM can be launched)
            changeMode(0);
        } else {
            var first  = text.indexOf("\""          );
            var second = text.indexOf("\"", first +1);
            var third  = text.indexOf("\"", second+1);
            var fourth = text.indexOf("\"", third +1);
            var curPhase = text.substring(first+1, second);
            var projName = text.substring(third+1, fourth);

            browserToServer("timLoad|capsMode|"+projName+":"+curPhase+"#<stealLock>|");
        }

        return;
        
    // check to see if another error occurred during loading
    } else if (tokens[1].substring(0,7) == "ERROR::") {
        alert(tokens[1]);

        // return the done menu to its original state
        var button = document.getElementById("doneMenuBtn");
        button["innerHTML"] = "Done";
        button.style.backgroundColor = null;

        var menu = document.getElementById("myDoneMenu");
        while (menu.firstChild != null) {
            menu.removeChild(menu.firstChild);
        }

        // return the solve button to its original state
        button = document.getElementById("solveButton");
        button["innerHTML"] = "Up to date";
        button.style.backgroundColor = null;

        // return to "main" mode (which is the only mode from which
        //    this TIM can be launched)
        changeMode(0);

    // otherwise get the intent phrase
    } else {
        wv.capsProj   = tokens[1];
        wv.capsPhase  = tokens[2];

        if (tokens[3] != "<haveIntent>") {
            wv.capsIntent = tokens[3];

            wv.capsIntent = prompt("Working on Phase "+wv.capsPhase+" of Project \""+wv.capsProj+"\"\n"+
                                   "What is the intent of the Phase?");
            if (wv.capsIntent === null) {
                wv.capsIntent = "<none>";
            }
            browserToServer("timMesg|capsMode|intent|"+wv.capsIntent+"|");
        }

        // update the legend on the Solve button
        button = document.getElementById("solveButton");
        button["innerHTML"] = wv.capsProj+":"+wv.capsPhase;

        // get the current Cvals so that the TreeWindow gets updated
        browserToServer("getPmtrs|");
        browserToServer("timMesg|capsMode|getCvals|13|");

        rebuildTreeWindow(97);

        // get the updated list of files
        browserToServer("getFilenames|");
    }
};


//
// callback when timSave returns
//
capsMode.timSaveCB = function (text) {
    alert("in capsMode.timSaveCB: "+text);

};


//
// callback when timQuit returns
//
capsMode.timQuitCB = function (text) {
    // alert("in capsMode.timQuitCB: "+text);

    // change solve button legend
    var button = document.getElementById("solveButton");
    button["innerHTML"] = "Up to date";
    button.style.backgroundColor = null;

    // restore the Done menu to its previous state
    var button = document.getElementById("doneMenuBtn");
    button["innerHTML"] = "Done";
    button.style.backgroundColor = null;

    var menu = document.getElementById("myDoneMenu");
    while (menu.firstChild != null) {
        menu.removeChild(menu.firstChild);
    }

    var button;

    button = document.createElement("input");
    button.type    = "button";
    button.title   = "Save and exit";
    button.value   = "Save";
    button.onclick = cmdSave;
    menu.appendChild(button);

    button = document.createElement("input");
    button.type    = "button";
    button.title   = "Quit and exit";
    button.value   = "Quit";
    button.onclick = cmdQuit;
    menu.appendChild(button);

    // toggle between hiding and showing the File menu contents
    document.getElementById("myDoneMenu").classList.toggle("showDoneMenu");

    // tell the tim to quit
    browserToServer("timQuit|capsMode|");

    // reset the server name and update the entries shown in the TreeWindow
    wv.server = "serveESP";
    rebuildTreeWindow(99);

    // change back to main mode
    changeMode(0);
};


//
// callback when timMesg returns
//
capsMode.timMesgCB = function (text) {
    alert("in capsMode.timMesgCB: "+text);

    var button = document.getElementById("solveButton");
    var tokens = text.split("|");

    // draw
    if (tokens[1] == "draw") {

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


//
// callback when "CAPS->UpdateIntent" is pressed
//
var cmdCapsUpdateIntent = function () {
    // alert("in cmdCapsUpdateIntent()");

    // close the Done menu
    document.getElementById("myDoneMenu").classList.remove("showDoneMenu");

    wv.capsIntent = prompt("Working on Phase "+wv.capsPhase+" of Project \""+wv.capsProj+"\"\n"+
                           "What is the intent of the Phase?", wv.capsIntent);
    if (wv.capsIntent === null) {
        wv.capsIntent = "<none>";
    }
    browserToServer("timMesg|capsMode|intent|"+wv.capsIntent+"|");
};


//
// callback when "CAPS->CommitPhase" is pressed
//
var cmdCapsCommitPhase = function () {
    // alert("in cmdCapsCommitPhase()");

    // close the current Phase
    browserToServer("timMesg|capsMode|commit|");

    postMessage("Committing \"" + wv.capsProj + ":" + wv.capsPhase+ "\"");
};


//
// callback when "CAPS->QuitPhase" is pressed
//
var cmdCapsQuitPhase = function () {
    // alert("in cmdCapsQuitPhase()");

    // close the current Phase (with remove)
    browserToServer("timMesg|capsMode|quit|");

    postMessage("Quitting \"" + wv.capsProj + ":" + wv.capsPhase+ "\"");

    // forget the Phase (since we really do not want to propose starting
    //    from a non-existant Phase
    wv.capsPhase = undefined;

    // remove the .csm file from the session (since it will be deleted as
    //    part of caps_close)
    alert("Quitting a Phase removes the .csm file from the session.\nUse File->Open to reload the .csm file.");

    browserToServer("new|");

    wv.filenames = capsMode.filenames;
    wv.nchanges  = 0;

    // remove from editor
    if (wv.codeMirror !== undefined) {
        wv.codeMirror.toTextArea();
        wv.codeMirror = undefined;
    }

    cval   = new Array();
    pmtr   = new Array();
    brch   = new Array();
    sgData = {};
};


//
// callback when "CAPS->SuspendPhase" is pressed
//
var cmdCapsSuspendPhase = function () {
    // alert("in cmdCapsSuspendPhase()");

    // close the current Phase (with remove)
    browserToServer("timMesg|capsMode|suspend|");

    postMessage("Suspending \"" + wv.capsProj + ":" + wv.capsPhase+ "\"");

    // forget the Phase (since we really do not want to propose starting
    //    from a non-existant Phase
    wv.capsPhase = undefined;
};


//
// callback when "CAPS->ListPhases" is pressed
//
var cmdCapsListPhases = function () {
    // alert("in cmdCapsListPhases()");

    // close the Done menu
    document.getElementById("myDoneMenu").classList.remove("showDoneMenu");

    // get the phases
    browserToServer("timMesg|capsMode|listPhases|");
};


//
// callback when "CAPS->ListHistory" is pressed
//
var cmdCapsListHistory = function () {
    // alert("in cmdCapsListHistory()");

    // close the Done menu
    document.getElementById("myDoneMenu").classList.remove("showDoneMenu");

    var myCval = prompt("What Caps Value do you want the history of?");
    if (myCval !== null) {

        // get the history
        browserToServer("timMesg|capsMode|listHistory|"+myCval+"|");
    }
};


//
// callback when "CAPS->ListAnalyses" is pressed
//
var cmdCapsListAnalyses = function () {
    // alert("in cmdCapsListAnalyses()");

    // close the Done menu
    document.getElementById("myDoneMenu").classList.remove("showDoneMenu");

    // get the analyses
    browserToServer("timMesg|capsMode|listAnalyses|");
};


//
// callback when "CAPS->ListBounds" is pressed
//
var cmdCapsListBounds = function () {
    // alert("in cmdCapsListBounds()");

    // close the Done menu
    document.getElementById("myDoneMenu").classList.remove("showDoneMenu");

    // get the bounds
    browserToServer("timMesg|capsMode|listBounds|");
};


//
// callback when "CAPS->ListAttachments" is pressed
//
var cmdCapsListAttachments = function () {
    alert("in cmdCapsListAttachments()");

    // close the Done menu
    document.getElementById("myDoneMenu").classList.remove("showDoneMenu");

    browserToServer("timMesg|capsMode|attachList|");
};


//
// callback when "CAPS->AttachFile" is pressed
//
var cmdCapsAttachFile = function () {
    alert("in cmdCapsAttachFile()");

    // close the Done menu
    document.getElementById("myDoneMenu").classList.remove("showDoneMenu");

    // attach a file
    var myAttach = prompt("Enter name of .txt, .png, .pdf, or .html file to attach");

    if (myAttach !== null) {
        browserToServer("timMesg|capsMode|attachFile|" + myAttach + "|");
    }
};


//
// callback when "CAPS->OpenAttachment" is pressed
//
var cmdCapsOpenAttachment = function () {
    alert("in cmdCapsOpenAttachment()");

    // close the Done menu
    document.getElementById("myDoneMenu").classList.remove("showDoneMenu");

    var myAttach = prompt("What attachment do you want to open?");

    if (myAttach !== null) {
        browserToServer("timMesg|capsMode|attachOpen|" + myAttach + "|");
    }
};
