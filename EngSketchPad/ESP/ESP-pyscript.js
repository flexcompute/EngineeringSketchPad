// ESP-pyscript.js implements python functions for the Engineering Sketch Pad (ESP)
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

// interface functions that ESP-pyscript.js provides
//    pyscript.launch()
//    pyscript.cmdUndo()
//    pyscript.cmdSolve()
//    pyscript.cmdSave()
//    pyscript.cmdQuit()
//
//    pyscript.cmdHome()                    not provided
//    pyscript.cmdLeft()                    not provided
//    pyscript.cmdRite()                    not provided
//    pyscript.cmdBotm()                    not provided
//    pyscript.cmdTop()                     not provided
//    pyscript.cmdIn()                      not provided
//    pyscript.cmdOut()                     not provided
//
//    pyscript.mouseDown(e)                 not provided
//    pyscript.mouseMove(e)                 not provided
//    pyscript.mouseUp(e)                   not provided
//    pyscript.mouseWheel(e)                not provided
//    pyscript.mouseLeftCanvas(e)           not provided
//
//    pyscript.keyPress(e)
//    pyscript.keyDown(e)                   not provided
//    pyscript.keyUp(e)                     not provided
//    pyscript.keyPressPart1(myKeyPress)    not provided
//    pyscript.keyPressPart2(picking, gprim)not provided
//    pyscript.updateKeyWindow()
//
//    pyscript.timLoadCB(text)
//    pyscript.timSaveCB(text)
//    pyscript.timQuitCB(text)
//    pyscript.timMesgCB(text)

"use strict";


//
// name of TIM object
//
pyscript.name = "pyscript";


//
// callback when Pyscript is launched
//
pyscript.launch = function (e) {
    // alert("in pyscript.launch(e="+e+")");

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

    // change mode
    changeMode(12);

    document.getElementById("toolMenuBtn").hidden = true;
    document.getElementById("doneMenuBtn").hidden = true;

    // load the tim and run the python script
    browserToServer("timLoad|pyscript|"+name+"|");

    // change solve button legend
    var button = document.getElementById("solveButton");
    button["innerHTML"] = "EditingPyscript";
    button.style.backgroundColor = "#FFFF3F";
};


//
// callback when "Pyscript->Undo" is pressed (called by ESP.html)
//
pyscript.cmdUndo = function () {
    alert("pyscript.cmdUndo() is not implemented");
};


//
// callback when "solveButton" is pressed
//
pyscript.cmdSolve = function () {
    // alert("in pyscript.cmdSolve()");

    var button = document.getElementById("solveButton");
    if (button["innerHTML"] == "EditingPyscript") {
        alert("Press <Cancel> or <Save and run> button on top of editor");
    } else if (button["innerHTML"] == "RunningPyscript") {
        alert("Pyscript execution cannot be interrupted")
    } else {
        alert("This should never happen (in pyscript.csmSolve)")
    }
};


//
// callback when "Pyscript->Save" is pressed (called by ESP.html)
//
pyscript.cmdSave = function () {
    // alert("in pyscript.cmdSave()");

    // if there were changed made in the editor, send the new file back to the server
    var fileChanged = 0;
    var newFile = wv.codeMirror.getDoc().getValue();
    if (pyscript.curFile != newFile) {
        postMessage("\""+pyscript.filename+"\" being saved and executed");
        fileChanged = 1;

        // because of an apparent limit on the size of text
        //    messages that can be sent from the browser to the
        //    server, we need to send the new file back in
        //    pieces and then reassemble on the server
        var maxMessageSize = 800;

        var ichar    = 0;
        var part     = newFile.substring(ichar, ichar+maxMessageSize);

        browserToServer("timMesg|pyscript|fileBeg|"+pyscript.filename+"|"+part);
        ichar += maxMessageSize;

        while (ichar < newFile.length) {
            part = newFile.substring(ichar, ichar+maxMessageSize);
            browserToServer("timMesg|pyscript|fileMid|"+part);
            ichar += maxMessageSize;
        }

        browserToServer("timMesg|pyscript|fileEnd|");
    } else {
        postMessage("\""+pyscript.filename+"\" being executed");
    }

    // change solve button legend
    var button = document.getElementById("solveButton");
    button["innerHTML"] = "RunningPyscript";

    // remove the editor
    if (wv.codeMirror !== undefined) {
        wv.codeMirror.toTextArea();
        wv.codeMirror = undefined;
        wv.fileindx   = undefined;
    }

    // restore buttons to their original states
    document.getElementById("toolMenuBtn").hidden = false;
    document.getElementById("editorCommentBtn").hidden = false;
    document.getElementById("editorIndentBtn").hidden = false;
    document.getElementById("editorHintBtn").hidden = false;
    document.getElementById("editorDebugBtn").hidden = false;

    // change the Cancel and Save buttons to their original states
    var cancelBtn = document.getElementById("editorCancelBtn");
    cancelBtn.value   = "Cancel";
    cancelBtn.onclick = editorCancel;

    var saveBtn = document.getElementById("editorSaveBtn");
    saveBtn.value   = "Save";
    saveBtn.onclick = editorOk;

    // restore the WebViewer
    document.getElementById("editorForm").hidden = true;
    document.getElementById("WebViewer").hidden = false;

    // run the python file
    if (fileChanged == 0) {
        browserToServer("timMesg|pyscript|execute|");
    }
};


//
// callback when "Pyscript->Quit" is pressed (called by ESP.html)
//
pyscript.cmdQuit = function () {
    // alert("in pyscript.cmdQuit()");

    // change solve button legend
    var button = document.getElementById("solveButton");
    button["innerHTML"] = "Up to date";
    button.style.backgroundColor = null;

    // remove the editor
    if (wv.codeMirror !== undefined) {
        wv.codeMirror.toTextArea();
        wv.codeMirror = undefined;
        wv.fileindx   = undefined;
    }

    // restore buttons to their original states
    document.getElementById("toolMenuBtn").hidden = false;
    document.getElementById("editorCommentBtn").hidden = false;
    document.getElementById("editorIndentBtn").hidden = false;
    document.getElementById("editorHintBtn").hidden = false;
    document.getElementById("editorDebugBtn").hidden = false;

    // change the Cancel and Save buttons to their original states
    var cancelBtn = document.getElementById("editorCancelBtn");
    cancelBtn.value   = "Cancel";
    cancelBtn.onclick = editorCancel;

    var saveBtn = document.getElementById("editorSaveBtn");
    saveBtn.value   = "Save";
    saveBtn.onclick = editorOk;

    // restore the WebViewer
    document.getElementById("editorForm").hidden = true;
    document.getElementById("WebViewer").hidden = false;

    // execute quit in the tim
    browserToServer("timQuit|pyscript|0|");
};


//
// callback when "homeButton" is pressed (calles by ESP.html)
//
//pyscript.cmdHome = function () {
//    main.cmdHome();
//};


//
// callback when "leftButton" is pressed (calles by ESP.html)
//
//pyscript.cmdLeft = function () {
//    main.cmdLeft();
//};


//
// callback when "riteButton" is pressed (calles by ESP.html)
//
//pyscript.cmdRite = function () {
//    main.cmdRite();
//};


//
// callback when "botmButton" is pressed (calles by ESP.html)
//
//pyscript.cmdBotm = function () {
//    main.cmdBotm();
//};


//
// callback when "topButton" is pressed (calles by ESP.html)
//
//pyscript.cmdTop = function () {
//    main.cmdTop();
//};


//
// callback when "inButton" is pressed (calles by ESP.html)
//
//pyscript.cmdIn = function () {
//    main.cmdIn();
//};


//
// callback when "outButton" is pressed (calles by ESP.html)
//
//pyscript.cmdOut = function () {
//    main.cmdOut();
//};


//
// callback when any mouse is pressed in canvas
//
//pyscript.mouseDown = function(e) {
//    main.mouseDown(e);
//};


//
// callback when any mouse moves in canvas
//
//pyscript.mouseMove = function(e) {
//    main.mouseMove(e);
//};


//
// callback when the mouse is released in canvas
//
//pyscript.mouseUp = function(e) {
//    main.mouseUp(e);
//};


//
// callback when the mouse wheel is rolled in canvas
//
//pyscript.mouseWheel = function(e) {
//    main.mouseWheel(e);
//};


//
// callback when the mouse leaves the canvas
//
//pyscript.mouseLeftCanvas = function(e) {
//    main.mouseLeftCanvas(e);
//};


//
// callback when a key is pressed
//
pyscript.keyPress = function (e) {
    // do nothing
};


//
// callback when an arror... or shift is pressed (needed for Chrome)
//
//pyscript.keyDown = function (e) {
//    main.keyDown(e);
//};


//
// callback when a shift is released (needed for Chrome)
//
//pyscript.keyUp = function (e) {
//    main.keyUp(e);
//};


//
// callback for first part of a keypress that is not recognized by wvUpdateUI
//
//pyscript.keyPressPart1 = function(myKeyPress) {
//    // alert("in pyscript.keyPressPart1(myKeyPress="+myKeyPress+")");
//    return 0;
//};


//
// callback for second part of a keypress that is not recognized by wvUpdateUI
//
//pyscript.keyPressPart2 = function(picking, gprim) {
//    // alert("in pyscript.keyPressPart2(picking="+picking+"   gprim="+gprim+")");
//};


//
// function to update the key window
//
//pyscript.updateKeyWindow = function () {
//};


//
// callback when timLoad returns
//
pyscript.timLoadCB = function (text) {
    // alert("in pyscript.timLoadCB: "+text.substring(7));

    if (text.substring(7,14) == "ERROR::") {
        alert(text.substring(7));

        pyscript.cmdQuit();

    } else {
        var tokens = text.split("|");
        if (tokens[2] == "") {
            return;
        }

        // remove former editor
        if (wv.codeMirror !== undefined && wv.myRole != 0) {
            wv.codeMirror.toTextArea();
            wv.codeMirror = undefined;
            wv.fileindx   = undefined;
        }

        pyscript.filename = tokens[1];
        pyscript.curFile  = tokens[2];
        for (var i = 3; i < tokens.length; i++) {
            pyscript.curFile += "|" + tokens[i];
        }

        var editorFilename = document.getElementById("editorFilename");
        editorFilename["innerHTML"] = "Contents of: "+pyscript.filename;

        var editorTextArea = document.getElementById("editorTextArea");

        editorTextArea.cols  = 84;
        editorTextArea.rows  = 25;
        editorTextArea.value = pyscript.curFile;

        // hide the WebViewer
        var webViewer     = document.getElementById("WebViewer");
        webViewer.hidden  = true;

        // unhide so that CodeMirror initialization will work
        var editorForm    = document.getElementById("editorForm");
        editorForm.hidden = false;

        // initialize CodeMirror
        if (wv.codeMirror === undefined) {
            var options = {
                lineNumbers:                     true,
                mode:                        "python",
                smartIndent:                     true,
                indentUnit:                         3,
                scrollbarStyle:              "native",
                matchBrackets:                   true,
                styleActiveLine:     {nonEmpty: true},
                styleActiveSelected: {nonEmpty: true},
                theme:                       "simple"
            };

            wv.codeMirror = CodeMirror.fromTextArea(editorTextArea, options);

            // choose between editable and readonly
            if (wv.myRole == 0) {
                wv.codeMirror.setOption("readOnly", false);
            } else {
                wv.codeMirror.setOption("readOnly",  true);
            }

            // send changes made in the editor to other collaborative users
            CodeMirror.on(wv.codeMirror, "change",
                          function (instance, obj) {
                              if (wv.myRole == 0 && wv.numUsers > 1) {
                                  if (obj.removed == "") {
                                      if (obj.text.length == 2) {
                                          browserToServer("editor|add|"+obj.from.line+"|"+obj.from.ch+"|<cr>|");
                                      } else {
                                          browserToServer("editor|add|"+obj.from.line+"|"+obj.from.ch+"|"+obj.text+"|");
                                      }
                                  } else {
                                      browserToServer("editor|del|"+obj.from.line+"|"+obj.from.ch+"|"+obj.to.line+"|"+obj.to.ch+"|");
                                  }
                              }
                          });

            var topRite = document.getElementById("trframe");
            var height  = Number(topRite.style.height.replace("px",""));
            wv.codeMirror.setSize(null, height-65);  // good value for nominal browser fonts
        }

        // clear any previous undo history and clipboard
        wv.codeMirror.clearHistory();
        wv.clipboard = "";

        // hide the Comment, Indent, Hint, and Debug  buttons
        document.getElementById("editorCommentBtn").hidden = true;
        document.getElementById("editorIndentBtn").hidden = true;
        document.getElementById("editorHintBtn").hidden = true;
        document.getElementById("editorDebugBtn").hidden = true;

        // change the Cancel and Save buttons
        var cancelBtn = document.getElementById("editorCancelBtn");
        cancelBtn.value   = "Cancel";
        cancelBtn.onclick = pyscript.cmdQuit;

        var saveBtn = document.getElementById("editorSaveBtn");
        saveBtn.value   = "Save and run";
        saveBtn.onclick = pyscript.cmdSave;
    }
};


//
// callback when timSave returns
//
pyscript.timSaveCB = function (text) {
    // alert("in pyscript.timSaveCB: "+text);
};


//
// callback when timQuit returns
//
pyscript.timQuitCB = function (text) {
    // alert("in pyscript.timQuitCB: "+text);

    document.getElementById("toolMenuBtn").hidden = false;
    document.getElementById("doneMenuBtn").hidden = true;

    // change solve button legend
    var button = document.getElementById("solveButton");
    button["innerHTML"] = "Up to date";
    button.style.backgroundColor = null;

    changeMode(0);
};


//
// callback when timMesg returns
//
pyscript.timMesgCB = function (text) {
    // alert("in pyscript.timMesgCB(text="+text+")")

    var tokens = text.split("|");

    // fileEnd
    if        (tokens[1] == "fileEnd") {
        browserToServer("timMesg|pyscript|execute|")

    // execute
    } else if (tokens[1] == "execute") {
//$$$        alert("in pyscript.timMesgCB: wv.overlay="+wv.overlay+"   wv.curTool.name="+wv.curTool.name);

        // update the display
        browserToServer("timDraw|pyscript|");
        browserToServer("getBrchs|");
        browserToServer("getPmtrs|");

        // get stderr if it exists
        browserToServer("timMesg|pyscript|stderr|");

        // quit
        browserToServer("timQuit|pyscript|1|");

        // change solve button legend
        var button = document.getElementById("solveButton");
        button["innerHTML"] = "Up to date";
        button.style.backgroundColor = null;

    // stdout
    } else if (tokens[1] == "stdout") {
        if (tokens[2].length > 0) {
            postMessage( "Output of pyscipt:\n"
                        +"******************\n"
                        + text.substring(14)
                        +"******************");
        } else {
            postMessage("Pyscript program created no output");
        }

    // stderr
    } else if (tokens[1] == "stderr") {
        if (tokens[2].length > 0) {
            alert( "Error generated by pyscript program:\n"
                   + text.substring(14, text.length-1));
            postMessage("\""+pyscript.filename+"\" completed with errors");
        } else {
            postMessage("\""+pyscript.filename+"\" completed successfully");
        }

        changeMode(0);

    // quit
    } else if (tokens[1] == "quit") {
        cmdQuit();

    // unknown message type
//    } else {
//        alert("unknown message type: "+ text);
    }
};

// /////////////////////////////////////////////////////////////////////

