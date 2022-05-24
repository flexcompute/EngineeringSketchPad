// ESP-sketch.js implements Sketcher functions for the Engineering Sketch Pad (ESP)
// written by John Dannenhoffer and Bob Haimes

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

// interface functions that ESP-sketch.js provides
//    sketch.launch()
//    sketch.cmdUndo()
//    sketch.cmdSolve()
//    sketch.cmdSave()
//    sketch.cmdQuit()
//
//    sketch.cmdHome()
//    sketch.cmdLeft()
//    sketch.cmdRite()
//    sketch.cmdBotm()
//    sketch.cmdTop()
//    sketch.cmdIn()
//    sketch.cmdOut()
//
//    sketch.mouseDown(e)
//    sketch.mouseMove(e)
//    sketch.mouseUp(e)
//
//    sketch.keyPress(e)

// functions associated with the Sketcher
//    sketcherInitialize()
//    sketcherSolvePre()
//    sketcherSolvePost(text)
//    sketcherSave()
//    sketcherUndo()
//    sketcherSaveUndo()
//    sketcherLoad()
//
//    sketcherUpdateSegData(iseg)
//    sketcherDraw()
//    sketcherDrawBezier(context, Xbezier, Ybezier)
//    sketcherDrawSpline(context, Xspline, Yspline)
//    sketcherGetClosestPoint()
//    sketcherGetClosestSegment()
//    sketcherGetClosestConstraint()

//    sketch.ibrch                     // Branch index that launched Sketcher
//    sketch.mode                      // 0 initiaizing
//                                     // 1 drawing
//                                     // 2 setting curvature
//                                     // 3 constraining
//                                     // 4 setting width
//                                     // 5 setting depth
//                                     // 6 solved
//    sketch.movingPoint = -1;         // moving point while constraining (or -1 for none)
//    sketch.basePoint = undefined;    // base Point for W or D constraint
//    sketch.skbegX    = undefined;    // expression for X from SKBEG
//    sketch.skbegY    = undefined;    // expression for Y from SKBEG
//    sketch.skbegZ    = undefined;    // expression for Z from SKBEG
//    sketch.relative  = 0;            // =1 if relative   from SKBEG
//    sketch.suggest   = undefined;    // suggested changes
//    sketch.heads[]                   // last  Point in a chain of Segments
//    sketch.tails[]                   // first Point in a chain of Segments
//    sketch.halo      =  5;           // pixels for determining "closeness" while drawing
//    sketch.offTop    =  0;           // offset to upper-left corner of the canvas
//    sketch.offLeft   =  0;

//    sketch.scale     = undefined;    // drawing scale factor
//    sketch.xorig     = undefined;    // x-origin for drawing
//    sketch.yorig     = undefined;    // y-origin for drawing

//    sketch.pnt.x[   i]               // x-coordinate
//    sketch.pnt.y[   i]               // y-coordinate
//    sketch.pnt.lbl[ i]               // constraint label

//    sketch.seg.type[i]               // Segment type (L, C, S, or B)
//    sketch.seg.ibeg[i]               // Point at beginning
//    sketch.seg.iend[i]               // Point at end
//    sketch.seg.xm[  i]               // x-coordinate at middle
//    sketch.seg.ym[  i]               // y-coordinate at middle
//    sketch.seg.lbl[ i]               // constraint label
//
//    sketch.seg.dip[ i]               // dip (>0 for convex, <0 for concave)
//    sketch.seg.xc[  i]               // x-coordinate at center
//    sketch.seg.yc[  i]               // y-coordinate at center
//    sketch.seg.rad[ i]               // radius (screen coordinates)
//    sketch.seg.tbeg[i]               // theta at beginning (radians)
//    sketch.seg.tend[i]               // theta at end       (radians)

//    sketch.var.name[ i]              // variable name
//    sketch.var.value[i]              // variable value

//    sketch.con.type[  i]             // constraint type (X, Y, P, T, A for Point)
//                                   //                 (H, V, I, L, R, or S for Segment)
//    sketch.con.index1[i]             // primary Point or Segment index
//    sketch.con.index2[i]             // secondary Point or Segment index
//    sketch.con.value[ i]             // (optional) constraint value

//    sketch.maxundo                   // maximum   undos in Sketcher
//    sketch.nundo                     // number of undos in Sketcher
//    sketch.undo[i].mode              // undo snapshot (or undefined)
//    sketch.undo[i].pnt
//    sketch.undo[i].seg
//    sketch.undo[i].var
//    sketch.undo[i].con

"use strict";


//
// name of TIM object
//
sketch.name = "sketch";


//
// callback when "Enter Sketcher" is pressed in editBrchForm (called by ESP.html)
//
sketch.launch = function () {
    // alert("sketch.launch()");

    // close the Tools menu
    var menu = document.getElementsByClassName("toolMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showToolMenu")) {
            menu[i].classList.remove(  "showToolMenu");
        }
    }

    sketch.ibrch = -1;
    var name     = "";

    // get the Tree node
    try {
        var id    = wv.menuEvent["target"].id;
        var inode = Number(id.substring(4,id.length-4));

        // get the Branch name
        name = myTree.name[inode].replace(/\u00a0/g, "").replace(/>/g, "");

        // get the Branch index
        for (var jbrch = 0; jbrch < brch.length; jbrch++) {
            if (brch[jbrch].name == name) {
                sketch.ibrch = jbrch + 1;
            }
        }
        if (sketch.ibrch <= 0) {
            alert(name+" not found");
            return;
        }

    /* if there is no Tree node, this implies that the sketecher
       was launched from the menu */
    } catch (err) {
        /* make a list of the Branches that contain a PATVAR */
        var nlist  = 0;
        var list   = [];
        var iskbeg = -1;
        for (var ibrch = 0; ibrch < brch.length; ibrch++) {
            if        (brch[ibrch].type == "skbeg") {
                iskbeg = ibrch;
            } else if (brch[ibrch].type == "skvar") {
                list[nlist] = iskbeg;
                nlist++;
            }
        }

        /* if there are PATBEG statements, find out from the user 
           which one should be editted */
        if (nlist > 0) {
            var mesg = "Enter\n0 to create new Sketch at end\n";

            for (var ilist = 0; ilist < nlist; ilist++) {
                mesg += (ilist+1) + " to edit Sketch starting at "+brch[list[ilist]].name+"\n";
            }

            var ans = prompt(mesg, "0");
            if (ans === null) {
                return;
            } else if (ans == 0) {
                sketch.ibrch = brch.length + 1;
                postMessage("Adding new sketch at (0,0,0) to end of Branches");
                browserToServer("newBrch|"+(sketch.ibrch-1)+"|skbeg|0|0|0|1||||||");
                browserToServer("getBrchs|");
                browserToServer("loadSketch|"+sketch.ibrch+"|");
                wv.nchanges++;
                changeMode(-1);
                return;
            } else if (ans > 0 && ans <= nlist) {
                sketch.ibrch = list[ans-1] + 1;
            } else {
                alert("Answer "+ans+" must be between 0 and "+nlist);
                return;
            }

        /* otherwise, put a PATBEG/PATEND after the last Branch */
        } else {
            sketch.ibrch = brch.length + 1;
            postMessage("Adding new sketch at (0,0,0) to end of Branches");
            browserToServer("newBrch|"+(sketch.ibrch-1)+"|skbeg|0|0|0|1||||||");
            browserToServer("getBrchs|");
            browserToServer("loadSketch|"+sketch.ibrch+"|");
            wv.nchanges++;
            changeMode(-1);
            return;
        }
    }
    
    // unhighlight the first column of the Tree
    unhighlightColumn1();

    // send the message to the server
    browserToServer("loadSketch|"+sketch.ibrch+"|");

    // inactivate buttons until build is done
    changeMode(-1);

    // initialize the Sketcher
    sketcherInitialize();
};


//
// callback when "undoButton" is pressed
//
sketch.cmdUndo = function () {
    // alert("in cmdUndo()");

    if (wv.server != "serveCSM") {
        alert("cmdUndo is not implemented for "+wv.server);
    } else {
        sketcherUndo();
    }
};


//
// callback when "solveButton" is pressed
//
sketch.cmdSolve = function () {
    // alert("in cmdSolve()");

    var button  = document.getElementById("solveButton");
    var buttext = button["innerHTML"];

    if (buttext == "Press to Solve") {
        if (wv.curMode != 8) {
            alert("Command disabled,  Press 'Cancel' or 'OK' first");
            return;
        }

        // save undo info
        sketcherSaveUndo();

        // send solve request to the server
        sketcherSolvePre();

    } else if (buttext == "Initializing...") {
        alert("The Sketcher is initializing.  Please be patient");

    } else if (buttext == "Drawing...") {
        alert("Use \"L\" to add a line segment\n" +
              "    \"C\" to add a circular arc\n" +
              "    \"S\" to add a spline point\n" +
              "    \"B\" to add a Bezier point\n" +
              "    \"Z\" to add a zero-length segment, or\n" +
              "    \"O\" to leave the Sketch open");

    } else if (buttext == "Setting R...") {
        alert("Press left mouse to set radius");

    } else if (buttext == "Constraining...") {

        // send solve request to the server
        sketcherSolvePre();

    } else if (buttext == "nothing") {
        alert("If at a point, use\n" +
              "    \"X\" to set the X-coordinate\n" +
              "    \"Y\" to set the Y-coordinate\n" +
              "    \"P\" to set perpendicularity\n" +
              "    \"T\" to set tangency\n" +
              "    \"A\" to set angle (positive to left)\n" +
              "    \"W\" to set width (dx) to another point\n" +
              "    \"D\" to set depth (dy) to another point\n" +
              "    \"<\" to delete a constraint\n" +
              "If at a cirarc segment, use\n" +
              "    \"R\" to set radius\n" +
              "    \"S\" to set sweep ccw sweep angle\n" +
              "    \"X\" to set X-coordinate of center\n" +
              "    \"Y\" to set Y-coordinate of center\n" +
              "    \"<\" to delete\n" +
              "If at a line segment, use\n" +
              "    \"H\" to set a horizontal constraint\n" +
              "    \"V\" to set a vertical constraint\n", +
              "    \"I\" to set inclination (from right horizontal)\n" +
              "    \"L\" to set length\n"+
              "    \"<\" to delete\n" +
              "If at a W or D constraint, use\n" +
              "    \"<\" to delete");

    } else if (buttext == "Setting W...") {
        alert("Move to another point and press left mouse to set width (dx)");

    } else if (buttext == "Setting D...") {
        alert("move to another point and press left mouse to set depth (dy)");

    } else {
        alert("Unexpected button text:\""+buttext+"\"");
    }
};


//
// callback when "Sketch->Save" is pressed (called by ESP.html)
//
sketch.cmdSave = function () {
    // alert("in sketch.cmdSave()");

    // close the Sketch menu
    var menu = document.getElementsByClassName("doneMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showDoneMenu")) {
            menu[i].classList.remove(  "showDoneMenu");
        }
    }

    if (wv.server != "serveCSM") {
        alert("cmdSave is not implemented for "+wv.server);
        return;

    } else if (wv.curMode == 8) {
        sketcherSave();
        wv.brchStat = 0;

    } else {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;

    }
};


//
// callback when "Sketch->Quit" is pressed (called by ESP.html)
//
sketch.cmdQuit = function () {
    // alert("in sketch.cmdQuit()");

    // close the Sketch menu
    var menu = document.getElementsByClassName("doneMenu-contents");
    for (var i = 0; i < menu.length; i++) {
        var openMenu = menu[i];
        if (menu[i].classList.contains("showDoneMenu")) {
            menu[i].classList.remove(  "showDoneMenu");
        }
    }

    if (wv.curMode != 8) {
        alert("Command disabled.  Press 'Cancel' or 'OK' first");
        return;
    } else if (wv.server != "serveCSM") {
        alert("cmdQuit is not implemented for "+wv.server);
        return;
    }

    changeMode(0);
    activateBuildButton();
};


//
// callback when "homeButton" is pressed
//
sketch.cmdHome = function () {
    var npnt = sketch.pnt.x.length;
    var nseg = sketch.seg.type.length;

    if (npnt <= 1) {
        alert("Not enough data to rescale");
        return;
    } else if (sketch.scale === undefined || isNaN(sketch.scale)) {
        alert("At least one length or radius must be set before rescaling  (sketch.scale="+sketch.scale+")");
        return;
    } else if (sketch.xorig === undefined || isNaN(sketch.xorig)) {
        alert("At lest one \"X\" must be set before rescaling  (sketch.xorig="+sketch.xorig+")");
        return;
    } else if (sketch.yorig === undefined || isNaN(sketch.yorig)) {
        alert("At lest one \"Y\" must be set before rescaling  (sketch.yorig="+sketch.yorig+")");
        return;
    }

    // get the extrema of the data
    var xmin = +1.0e+10;
    var xmax = -1.0e+10;
    var ymin = +1.0e+10;
    var ymax = -1.0e+10;

    for (var ipnt = 0; ipnt < npnt; ipnt++) {
        xmin = Math.min(xmin, sketch.scale * (sketch.pnt.x[ipnt] - sketch.xorig));
        xmax = Math.max(xmax, sketch.scale * (sketch.pnt.x[ipnt] - sketch.xorig));
        ymin = Math.min(ymin, sketch.scale * (sketch.yorig - sketch.pnt.y[ipnt]));
        ymax = Math.max(ymax, sketch.scale * (sketch.yorig - sketch.pnt.y[ipnt]));
    }

    // get the size of the canvas
    var canvas = document.getElementById("sketcher");
    var width  = canvas.clientWidth;
    var height = canvas.clientHeight;

    var test1 = (xmax - xmin) / width;
    var test2 = (ymax - ymin) / height;

    // set up sizing so that Sketch fills middle 50 percent of the window
    var newScale = 2 * Math.max(test1, test2);
    var newXorig = (width  - (xmin + xmax) / newScale) / 2;
    var newYorig = (height + (ymin + ymax) / newScale) / 2;

    // convert the Points by the new scale factors
    for (ipnt = 0; ipnt < npnt; ipnt++) {
        var xx = sketch.scale * (sketch.pnt.x[ipnt] - sketch.xorig);
        var yy = sketch.scale * (sketch.yorig - sketch.pnt.y[ipnt]);

//        sketch.pnt.x[ipnt] = Math.floor(newXorig + xx / newScale);
//        sketch.pnt.y[ipnt] = Math.floor(newYorig - yy / newScale);
        sketch.pnt.x[ipnt] = newXorig + xx / newScale;
        sketch.pnt.y[ipnt] = newYorig - yy / newScale;
    }

    for (var iseg = 0; iseg < nseg; iseg++) {
        sketch.seg.dip[iseg] = sketch.seg.dip[iseg] * sketch.scale / newScale;

        sketcherUpdateSegData(iseg);
    }

    sketch.scale = newScale;
    sketch.xorig = newXorig;
    sketch.yorig = newYorig;

    sketcherDraw();
};


//
// callback when "leftButton" is pressed
//
sketch.cmdLeft = function () {
    var npnt = sketch.pnt.x.length;
    var nseg = sketch.seg.type.length;

    if (npnt <= 1) {
        alert("Not enough data to rescale");
        return;
    } else if (sketch.scale === undefined || isNaN(sketch.scale)) {
        alert("At least one length or radius must be set before rescaling  (sketch.scale="+sketch.scale+")");
        return;
    } else if (sketch.xorig === undefined || isNaN(sketch.xorig)) {
        alert("At lest one \"X\" must be set before rescaling  (sketch.xorig="+sketch.xorig+")");
        return;
    } else if (sketch.yorig === undefined || isNaN(sketch.yorig)) {
        alert("At lest one \"Y\" must be set before rescaling  (sketch.yorig="+sketch.yorig+")");
        return;
    }

    var newXorig = sketch.xorig - 100;

    for (var ipnt = 0; ipnt < npnt; ipnt++) {
        var xx = sketch.scale * (sketch.pnt.x[ipnt] - sketch.xorig);
//        sketch.pnt.x[ipnt] = Math.floor(newXorig + xx / sketch.scale);
        sketch.pnt.x[ipnt] = newXorig + xx / sketch.scale;
    }

    for (var iseg = 0; iseg < nseg; iseg++) {
        sketcherUpdateSegData(iseg);
    }

    sketch.xorig = newXorig;

    sketcherDraw();
};


//
// callback when "riteButton" is pressed
//
sketch.cmdRite = function () {
    var npnt = sketch.pnt.x.length;
    var nseg = sketch.seg.type.length;

    if (npnt <= 1) {
        alert("Not enough data to rescale");
        return;
    } else if (sketch.scale === undefined || isNaN(sketch.scale)) {
        alert("At least one length or radius must be set before rescaling  (sketch.scale="+sketch.scale+")");
        return;
    } else if (sketch.xorig === undefined || isNaN(sketch.xorig)) {
        alert("At lest one \"X\" must be set before rescaling  (sketch.xorig="+sketch.xorig+")");
        return;
    } else if (sketch.yorig === undefined || isNaN(sketch.yorig)) {
        alert("At lest one \"Y\" must be set before rescaling  (sketch.yorig="+sketch.yorig+")");
        return;
    }

    var newXorig = sketch.xorig + 100;

    for (var ipnt = 0; ipnt < npnt; ipnt++) {
        var xx = sketch.scale * (sketch.pnt.x[ipnt] - sketch.xorig);
//        sketch.pnt.x[ipnt] = Math.floor(newXorig + xx / sketch.scale);
        sketch.pnt.x[ipnt] = newXorig + xx / sketch.scale;
    }

    for (var iseg = 0; iseg < nseg; iseg++) {
        sketcherUpdateSegData(iseg);
    }

    sketch.xorig = newXorig;

    sketcherDraw();
};


//
// callback when "botmButton" is pressed
//
sketch.cmdBotm = function () {
    var npnt = sketch.pnt.x.length;
    var nseg = sketch.seg.type.length;

    if (npnt <= 1) {
        alert("Not enough data to rescale");
        return;
    } else if (sketch.scale === undefined || isNaN(sketch.scale)) {
        alert("At least one length or radius must be set before rescaling  (sketch.scale="+sketch.scale+")");
        return;
    } else if (sketch.xorig === undefined || isNaN(sketch.xorig)) {
        alert("At lest one \"X\" must be set before rescaling  (sketch.xorig="+sketch.xorig+")");
        return;
    } else if (sketch.yorig === undefined || isNaN(sketch.yorig)) {
        alert("At lest one \"Y\" must be set before rescaling  (sketch.yorig="+sketch.yorig+")");
        return;
    }

    var newYorig = sketch.yorig + 100;

    for (var ipnt = 0; ipnt < npnt; ipnt++) {
        var yy = sketch.scale * (sketch.yorig - sketch.pnt.y[ipnt]);
//        sketch.pnt.y[ipnt] = Math.floor(newYorig - yy / sketch.scale);
        sketch.pnt.y[ipnt] = newYorig - yy / sketch.scale;
    }

    for (var iseg = 0; iseg < nseg; iseg++) {
        sketcherUpdateSegData(iseg);
    }

    sketch.yorig = newYorig;

    sketcherDraw();
};


//
// callback when "topButton" is pressed
//
sketch.cmdTop = function () {
    var npnt = sketch.pnt.x.length;
    var nseg = sketch.seg.type.length;

    if (npnt <= 1) {
        alert("Not enough data to rescale");
        return;
    } else if (sketch.scale === undefined || isNaN(sketch.scale)) {
        alert("At least one length or radius must be set before rescaling  (sketch.scale="+sketch.scale+")");
        return;
    } else if (sketch.xorig === undefined || isNaN(sketch.xorig)) {
        alert("At lest one \"X\" must be set before rescaling  (sketch.xorig="+sketch.xorig+")");
        return;
    } else if (sketch.yorig === undefined || isNaN(sketch.yorig)) {
        alert("At lest one \"Y\" must be set before rescaling  (sketch.yorig="+sketch.yorig+")");
        return;
    }

    var newYorig = sketch.yorig - 100;

    for (var ipnt = 0; ipnt < npnt; ipnt++) {
        var yy = sketch.scale * (sketch.yorig - sketch.pnt.y[ipnt]);
//        sketch.pnt.y[ipnt] = Math.floor(newYorig - yy / sketch.scale);
        sketch.pnt.y[ipnt] = newYorig - yy / sketch.scale;
    }

    for (var iseg = 0; iseg < nseg; iseg++) {
        sketcherUpdateSegData(iseg);
    }

    sketch.yorig = newYorig;

    sketcherDraw();
};


//
// callback when "inButton" is pressed
//
sketch.cmdIn = function () {
    var npnt = sketch.pnt.x.length;
    var nseg = sketch.seg.type.length;

    if (npnt <= 1) {
        alert("Not enough data to rescale");
        return;
    } else if (sketch.scale === undefined || isNaN(sketch.scale)) {
        alert("At least one length or radius must be set before rescaling  (sketch.scale="+sketch.scale+")");
        return;
    } else if (sketch.xorig === undefined || isNaN(sketch.xorig)) {
        alert("At lest one \"X\" must be set before rescaling  (sketch.xorig="+sketch.xorig+")");
        return;
    } else if (sketch.yorig === undefined || isNaN(sketch.yorig)) {
        alert("At lest one \"Y\" must be set before rescaling  (sketch.yorig="+sketch.yorig+")");
        return;
    }

    var newScale = 0.5 * sketch.scale;

    for (var ipnt = 0; ipnt < npnt; ipnt++) {
        var xx = sketch.scale * (sketch.pnt.x[ipnt] - sketch.xorig);
        var yy = sketch.scale * (sketch.yorig - sketch.pnt.y[ipnt]);

//        sketch.pnt.x[ipnt] = Math.floor(sketch.xorig + xx / newScale);
//        sketch.pnt.y[ipnt] = Math.floor(sketch.yorig - yy / newScale);
        sketch.pnt.x[ipnt] = sketch.xorig + xx / newScale;
        sketch.pnt.y[ipnt] = sketch.yorig - yy / newScale;
    }

    for (var iseg = 0; iseg < nseg; iseg++) {
        sketch.seg.dip[iseg] = sketch.seg.dip[iseg] * sketch.scale / newScale;

        sketcherUpdateSegData(iseg);
    }

    sketch.scale = newScale;

    sketcherDraw();
};


//
// callback when "outButton" is pressed
//
sketch.cmdOut = function () {
    var npnt = sketch.pnt.x.length;
    var nseg = sketch.seg.type.length;

    if (npnt <= 1) {
        alert("Not enough data to rescale");
        return;
    } else if (sketch.scale === undefined || isNaN(sketch.scale)) {
        alert("At least one length or radius must be set before rescaling  (sketch.scale="+sketch.scale+")");
        return;
    } else if (sketch.xorig === undefined || isNaN(sketch.xorig)) {
        alert("At lest one \"X\" must be set before rescaling  (sketch.xorig="+sketch.xorig+")");
        return;
    } else if (sketch.yorig === undefined || isNaN(sketch.yorig)) {
        alert("At lest one \"Y\" must be set before rescaling  (sketch.yorig="+sketch.yorig+")");
        return;
    }

    var newScale = 2.0 * sketch.scale;

    for (var ipnt = 0; ipnt < npnt; ipnt++) {
        var xx = sketch.scale * (sketch.pnt.x[ipnt] - sketch.xorig);
        var yy = sketch.scale * (sketch.yorig - sketch.pnt.y[ipnt]);

//        sketch.pnt.x[ipnt] = Math.floor(sketch.xorig + xx / newScale);
//         sketch.pnt.y[ipnt] = Math.floor(sketch.yorig - yy / newScale);
        sketch.pnt.x[ipnt] = sketch.xorig + xx / newScale;
        sketch.pnt.y[ipnt] = sketch.yorig - yy / newScale;
    }

    for (var iseg = 0; iseg < nseg; iseg++) {
        sketch.seg.dip[iseg] = sketch.seg.dip[iseg] * sketch.scale / newScale;

        sketcherUpdateSegData(iseg);
    }

    sketch.scale = newScale;

    sketcherDraw();
};


//
// callback when any mouse is pressed  (when wv.curMode==8)
//
sketch.mouseDown = function (e) {
    if (!e) var e = event;

    wv.startX   = e.clientX - sketch.offLeft - 1;
    wv.startY   = e.clientY - sketch.offTop  - 1;

    wv.button   = e.button;

                    wv.modifier  = 0;
    if (e.shiftKey) wv.modifier |= 1;
    if (e.altKey  ) wv.modifier |= 2;
    if (e.ctrlKey ) wv.modifier |= 4;

    // set default key presses that are bound to mouse clicks
    if        (sketch.mode == 1) {
        var evnt = new Object();
        evnt.keyCode  = 0;
        evnt.charCode = 108;   // "l"

        sketch.keyPress(evnt);
    } else if (sketch.mode == 2) {
        var evnt = new Object();
        evnt.keyCode  = 0;
        evnt.charCode = 108;   // "l"

        sketch.keyPress(evnt);
    } else if (sketch.mode == 3) {
        wv.cursorX = wv.startX;
        wv.cursorY = wv.startY;

        sketch.movingPoint = sketcherGetClosestPoint();
    } else if (sketch.mode == 4) {
        var evnt = new Object();
        evnt.keyCode  = 0;
        evnt.charCode = 124;   // "w"

        sketch.keyPress(evnt);
    } else if (sketch.mode == 5) {
        var evnt = new Object();
        evnt.keyCode  = 0;
        evnt.charCode = 100;   // "d"

        sketch.keyPress(evnt);
    }
};


//
// callback when the mouse moves
//
sketch.mouseMove = function (e) {
    if (!e) var e = event;

    wv.cursorX  = e.clientX - sketch.offLeft - 1;
    wv.cursorY  = e.clientY - sketch.offTop  - 1;

                    wv.modifier  = 0;
    if (e.shiftKey) wv.modifier |= 1;
    if (e.altKey  ) wv.modifier |= 2;
    if (e.ctrlKey ) wv.modifier |= 4;

    // if sketch.mode==1, then the "proposed" new Segment
    //    is taken care of in sketcherDraw
    if (sketch.mode == 1) {
        sketcherDraw();

    // if sketch.mode==3, then move the point at (wv.startX, wv.startY)
    //    to the current location
    } else if (sketch.mode == 3) {
        if (sketch.movingPoint >= 0) {
            sketch.pnt.x[sketch.movingPoint] = wv.cursorX;
            sketch.pnt.y[sketch.movingPoint] = wv.cursorY;

            for (var iseg = 0; iseg < sketch.seg.type.length; iseg++) {
                sketcherUpdateSegData(iseg);
            }

            sketcherDraw();
        }

    // if sketch.mode==4 or 5, then the "proposed" new W or D constraint
    //    is taken care of in sketcherDraw
    } else if (sketch.mode == 4 || sketch.mode == 5) {
        sketcherDraw();

    // if sketch.mode==2, then we need to use the cursor location to update
    //    the dip for the last Segment
    } else if (sketch.mode == 2) {
        var iseg = sketch.seg.type.length - 1;

        var xa = sketch.pnt.x[sketch.seg.ibeg[iseg]];
        var ya = sketch.pnt.y[sketch.seg.ibeg[iseg]];
        var xb = sketch.pnt.x[sketch.seg.iend[iseg]];
        var yb = sketch.pnt.y[sketch.seg.iend[iseg]];

        // put the cursor on the arc
        var xe  = wv.cursorX;
        var ye  = wv.cursorY;

        var D = (ya - ye) * (xb - xe) - (yb - ye) * (xa - xe);
        if (Math.abs(D) > 1e-6) {
            var s = ((xb - xa) * (xb - xe) - (ya - yb) * (yb - ye)) / D;

            var xc = (xa + xe + s * (ya - ye)) / 2;
            var yc = (ya + ye + s * (xe - xa)) / 2;

            var R = Math.sqrt((xc-xa) * (xc-xa) + (yc-ya) * (yc-ya));
            var L = Math.sqrt((xb-xa) * (xb-xa) + (yb-ya) * (yb-ya));

            if ((xe-xa)*(yb-ya) > (ye-ya)*(xb-xa)) {
                if ((xb-xa)*(yc-ya) > (yb-ya)*(xc-xa)) {
                    sketch.seg.dip[iseg] = -R + Math.sqrt(R*R - L*L/4);
                } else {
                    sketch.seg.dip[iseg] = -R - Math.sqrt(R*R - L*L/4);
                }
            } else {
                if ((xb-xa)*(yc-ya) > (yb-ya)*(xc-xa)) {
                    sketch.seg.dip[iseg] = +R + Math.sqrt(R*R - L*L/4);
                } else {
                    sketch.seg.dip[iseg] = +R - Math.sqrt(R*R - L*L/4);
                }
            }

            sketcherUpdateSegData(iseg);
        } else {
            sketch.seg.dip[ iseg] = 0;
            sketch.seg.xm[  iseg] = (xa + xb) / 2;
            sketch.seg.ym[  iseg] = (ya + yb) / 2;
            sketch.seg.xc[  iseg] = 0;
            sketch.seg.yc[  iseg] = 0;
            sketch.seg.rad[ iseg] = 0;
            sketch.seg.tbeg[iseg] = 0;
            sketch.seg.tend[iseg] = 0;
        }

        sketcherDraw();
    }
};


//
// callback when the mouse is released
//
sketch.mouseUp = function (e) {
    sketch.movingPoint = -1;
};


//
// process a key press in the Sketcher
//
sketch.keyPress = function (e) {
    // alert("sketch.keyPress");

    // record info about keypress
    wv.keyPress = e.charCode;
    wv.keyCode  = e.keyCode;

                    wv.modifier  = 0;
    if (e.shiftKey) wv.modifier |= 1;
    if (e.altKey  ) wv.modifier |= 2;
    if (e.ctrlKey ) wv.modifier |= 4;

    var canvas  = document.getElementById("sketcher");
    var context = canvas.getContext("2d");

    var npnt = sketch.pnt.x.length;
    var nseg = sketch.seg.type.length;
    var nvar = sketch.var.name.length;
    var ncon = sketch.con.type.length;

    var myKeyPress = String.fromCharCode(wv.keyPress);

    // 'ctrl-h' --- home
    if ((wv.keyPress == 104 && wv.modifier == 4 && wv.keyCode == 0) ||
        (wv.keyPress ==   8 && wv.modifier == 4 && wv.keyCode == 8)   ) {
        sketch.cmdHome();
        return;

    // '<Home>' -- home
    } else if (wv.keyPress == 0 && wv.keyCode == 36) {
        sketch.cmdHome();
        return;

    // 'ctrl-i' - zoom in (same as <PgUp> without shift)
    } else if ((wv.keyPress == 105 && wv.modifier == 4 && wv.keyCode ==  0) ||
               (wv.keyPress ==   9 && wv.modifier == 4 && wv.keyCode ==  9)   ) {
        sketch.cmdIn();
        return;

    // '<PgDn>' -- zoom in
    } else if (wv.keyPress == 0 && wv.keyCode == 33) {
        sketch.cmdIn();
        return;

    // '+' -- zoom in
    } else if (wv.keyPress == 43 && wv.modifier == 1) {
        sketch.cmdIn();
        return;

    // 'ctrl-o' - zoom out (same as <PgDn> without shift)
    } else if ((wv.keyPress == 111 && wv.modifier == 4 && wv.keyCode ==  0) ||
               (wv.keyPress ==  15 && wv.modifier == 4 && wv.keyCode == 15)   ) {
        sketch.cmdOut();
        return;

    // '<PgUp>' -- zoom out
    } else if (wv.keyPress == 0 && wv.keyCode == 34) {
        sketch.cmdOut();
        return;

    // '-' -- zoomout
    } else if (wv.keyPress == 45 && wv.modifier == 0) {
        sketch.cmdOut();
        return;

    // 'ctrl-r' - translate right
    } else if ((wv.keyPress == 114 && wv.modifier == 4 && wv.keyCode ==  0) ||
               (wv.keyPress ==  18 && wv.modifier == 4 && wv.keyCode == 18)   ) {
        sketch.cmdRite();
        return;

    // '<Right>' -- translate right
    } else if (wv.keyPress == 0 && wv.keyCode == 39) {
        sketch.cmdRite();
        return;

    // 'ctrl-l' - translate left
    } else if ((wv.keyPress == 108 && wv.modifier == 4 && wv.keyCode ==  0) ||
               (wv.keyPress ==  12 && wv.modifier == 4 && wv.keyCode == 12)   ) {
        sketch.cmdLeft();
        return;

    // '<Left>' -- translate left
    } else if (wv.keyPress == 0 && wv.keyCode == 37) {
        sketch.cmdLeft();
        return;

    // 'ctrl-t' - translate up
    } else if ((wv.keyPress == 116 && wv.modifier == 4 && wv.keyCode ==  0) ||
               (wv.keyPress ==  20 && wv.modifier == 4 && wv.keyCode == 20)   ) {
        sketch.cmdTop();
        return;

    // '<Up>' -- translate up
    } else if (wv.keyPress == 0 && wv.keyCode == 38) {
        sketch.cmdTop();
        return;

    // 'ctrl-b' - translate down
    } else if ((wv.keyPress ==  98 && wv.modifier == 4 && wv.keyCode ==  0) ||
               (wv.keyPress ==   2 && wv.modifier == 4 && wv.keyCode ==  2)   ) {
        sketch.cmdBotm();
        return;

    // '<Down>' -- translate down
    } else if (wv.keyPress == 0 && wv.keyCode == 40) {
        sketch.cmdBotm();
        return;

    // '@' - coordinates of cursor
    } else if (myKeyPress == "@" || myKeyPress == '2') {
        var xx = sketch.scale * (wv.cursorX - sketch.xorig);
        var yy = sketch.scale * (sketch.yorig - wv.cursorY);
        postMessage("Cursor is at ("+xx+","+yy+")");

    // '*' - rescale based upon best info
    } else if (myKeyPress == '*' || myKeyPress == '8') {
        var xold = sketch.scale * (wv.cursorX - sketch.xorig);
        var yold = sketch.scale * (sketch.yorig - wv.cursorY);

        var xnew = prompt("Enter approximate X location at cursor: ", xold);
        if (xnew !== null) {
            var ynew = prompt("Enter approximate Y location at cursor: ", yold);
            if (ynew !== null) {
                postMessage("Resetting scale based upon cursor location");

                var Lphys = Math.sqrt(xnew*xnew + ynew*ynew);
                var Lscr  = Math.sqrt(Math.pow(wv.cursorX-sketch.xorig,2)
                                     +Math.pow(wv.cursorY-sketch.yorig,2));
                sketch.scale = Lphys / Lscr;
            }
        }

    // sketch.mode 1 options ('l', 'c', 's', 'b', 'z', or 'o' to leave sketch open)
    } else if (sketch.mode == 1) {
        sketcherSaveUndo();

        // determine the coordinates (and constraints) for this placement
        var xnew = wv.cursorX;
        var ynew = wv.cursorY;
        var cnew = "";

        for (var ipnt = npnt-1; ipnt >= 0; ipnt--) {
            if (Math.abs(wv.cursorX-sketch.pnt.x[ipnt]) < sketch.halo) {
                xnew = sketch.pnt.x[ipnt];
                if (ipnt == npnt-1) {
                    cnew = "V";
                }
                break;
            }
        }
        for (ipnt = npnt-1; ipnt >= 0; ipnt--) {
            if (Math.abs(wv.cursorY-sketch.pnt.y[ipnt]) < sketch.halo) {
                ynew = sketch.pnt.y[ipnt];
                if (ipnt == npnt-1) {
                    cnew = "H";
                }
                break;
            }
        }

        // "l" --- add a line Segment
        if (myKeyPress == "l" || myKeyPress == "L") {
            var ibeg = npnt - 1;
            var iend = 0;

            // if this point is within halo of previous point, reject it
            if (Math.abs(sketch.pnt.x[npnt-1]-xnew) < sketch.halo &&
                Math.abs(sketch.pnt.y[npnt-1]-ynew) < sketch.halo   ) {
                postMessage("cannot create zero-length line segment\n use \"Z\" command instead");
                return;
            }

            // if this Point does not close the Sketch, add a new Point and its variables
            if (Math.abs(sketch.pnt.x[0]-xnew) > sketch.halo ||
                Math.abs(sketch.pnt.y[0]-ynew) > sketch.halo   ) {
                iend = npnt;
                sketch.pnt.x[  iend] = xnew;
                sketch.pnt.y[  iend] = ynew;
                sketch.pnt.lbl[iend] = "";

                sketch.var.name[ nvar  ] = sprintf("::x[%]", npnt+1);
                sketch.var.value[nvar  ] = sprintf("%",    xnew  );

                sketch.var.name[ nvar+1] = sprintf("::y[%]", npnt+1);
                sketch.var.value[nvar+1] = sprintf("%",    ynew  );
            }

            // add the new Segment
            sketch.seg.type[nseg] = "L";
            sketch.seg.ibeg[nseg] = ibeg;
            sketch.seg.iend[nseg] = iend;
            sketch.seg.xm[  nseg] = (sketch.pnt.x[ibeg] + sketch.pnt.x[iend]) / 2;
            sketch.seg.ym[  nseg] = (sketch.pnt.y[ibeg] + sketch.pnt.y[iend]) / 2;
            sketch.seg.lbl[ nseg] = "";

            sketch.seg.dip[ nseg] = 0;
            sketch.seg.xc[  nseg] = 0;
            sketch.seg.yc[  nseg] = 0;
            sketch.seg.rad[ nseg] = 0;
            sketch.seg.tbeg[nseg] = 0;
            sketch.seg.tend[nseg] = 0;

            // add a new "horizontal" or "vertical" constraint
            if        (cnew == "H") {
                sketch.con.type[  ncon] = "H";
                sketch.con.index1[ncon] = ibeg;
                sketch.con.index2[ncon] = iend;
                sketch.con.value[ ncon] = 0;

                sketch.seg.lbl[nseg] += "H";
            } else if (cnew == "V") {
                sketch.con.type[  ncon] = "V";
                sketch.con.index1[ncon] = ibeg;
                sketch.con.index2[ncon] = iend;
                sketch.con.value[ ncon] = 0;

                sketch.seg.lbl[nseg] += "V";
            }

            // if Sketch is closed, post a message and swith mode
            if (iend == 0) {
                sketch.mode = 3;
            }

        // "c" --- add a circular arc
        } else if (myKeyPress == "c" || myKeyPress == "C") {
            var ibeg = npnt - 1;
            var iend = 0;

            // if this point is within halo of previous point, reject it
            if (Math.abs(sketch.pnt.x[npnt-1]-xnew) < sketch.halo &&
                Math.abs(sketch.pnt.y[npnt-1]-ynew) < sketch.halo   ) {
                postMessage("cannot create zero-length circular arc");
                return;
            }

            // if this Point does not close the Sketch, add a new Point and its variables
            if (Math.abs(sketch.pnt.x[0]-xnew) > sketch.halo ||
                Math.abs(sketch.pnt.y[0]-ynew) > sketch.halo) {
                iend          = npnt;
                sketch.pnt.x[  iend] = xnew;
                sketch.pnt.y[  iend] = ynew;
                sketch.pnt.lbl[iend] = "";

                sketch.var.name[ nvar  ] = sprintf("::x[%]", npnt+1);
                sketch.var.value[nvar  ] = sprintf("%",    xnew  );

                sketch.var.name[ nvar+1] = sprintf("::y[%]", npnt+1);
                sketch.var.value[nvar+1] = sprintf("%",    ynew  );
            }

            // add the new Segment
            sketch.seg.type[nseg] = "C";
            sketch.seg.ibeg[nseg] = ibeg;
            sketch.seg.iend[nseg] = iend;
            sketch.seg.xm[  nseg] = (sketch.pnt.x[ibeg] + sketch.pnt.x[iend]) / 2;
            sketch.seg.ym[  nseg] = (sketch.pnt.y[ibeg] + sketch.pnt.y[iend]) / 2;
            sketch.seg.lbl[ nseg] = "";

            sketch.seg.dip[ nseg] = 0;
            sketch.seg.xc[  nseg] = 0;
            sketch.seg.yc[  nseg] = 0;
            sketch.seg.rad[ nseg] = 0;
            sketch.seg.tbeg[nseg] = 0;
            sketch.seg.tend[nseg] = 0;

            nvar = sketch.var.name.length;

            sketch.var.name[ nvar] = sprintf("::d[%]", sketch.seg.type.length);
            sketch.var.value[nvar] = "0";

            sketch.mode = 2;

        // "s" --- add spline Point
        } else if (myKeyPress == "s" || myKeyPress == "S") {
            var ibeg = npnt - 1;
            var iend = 0;

            // if this point is within halo of previous point, reject it
            if (Math.abs(sketch.pnt.x[npnt-1]-xnew) < sketch.halo &&
                Math.abs(sketch.pnt.y[npnt-1]-ynew) < sketch.halo   ) {
                postMessage("cannot create zero-length spline segment");
                return;
            }

            // if this Point does not close the Sketch, add a new Point and its variables
            if (Math.abs(sketch.pnt.x[0]-xnew) > sketch.halo ||
                Math.abs(sketch.pnt.y[0]-ynew) > sketch.halo   ) {
                iend = npnt;
                sketch.pnt.x[  iend] = xnew;
                sketch.pnt.y[  iend] = ynew;
                sketch.pnt.lbl[iend] = "";

                sketch.var.name[ nvar  ] = sprintf("::x[%]", npnt+1);
                sketch.var.value[nvar  ] = sprintf("%",    xnew  );

                sketch.var.name[ nvar+1] = sprintf("::y[%]", npnt+1);
                sketch.var.value[nvar+1] = sprintf("%",    ynew  );
            }

            // add the new Segment
            sketch.seg.type[nseg] = "S";
            sketch.seg.ibeg[nseg] = ibeg;
            sketch.seg.iend[nseg] = iend;
            sketch.seg.xm[  nseg] = (sketch.pnt.x[ibeg] + sketch.pnt.x[iend]) / 2;
            sketch.seg.ym[  nseg] = (sketch.pnt.y[ibeg] + sketch.pnt.y[iend]) / 2;
            sketch.seg.lbl[ nseg] = "";

            sketch.seg.dip[ nseg] = 0;
            sketch.seg.xc[  nseg] = 0;
            sketch.seg.yc[  nseg] = 0;
            sketch.seg.rad[ nseg] = 0;
            sketch.seg.tbeg[nseg] = 0;
            sketch.seg.tend[nseg] = 0;

            // if Sketch is closed, post a message and swith mode
            if (iend == 0) {
                sketch.mode = 3;
            }

        // "b" --- add Bezier Point
        } else if (myKeyPress == "b" || myKeyPress == "B") {
            var ibeg = npnt - 1;
            var iend = 0;

            // if this point is within halo of last point, reject it
            if (Math.abs(sketch.pnt.x[npnt-1]-xnew) < sketch.halo &&
                Math.abs(sketch.pnt.y[npnt-1]-ynew) < sketch.halo   ) {
                postMessage("cannot create zero-length bezier segment");
                return;
            }

            // if this Point does not close the Sketch, add a new Point and its variables
            if (Math.abs(sketch.pnt.x[0]-xnew) > sketch.halo ||
                Math.abs(sketch.pnt.y[0]-ynew) > sketch.halo   ) {
                iend = npnt;
                sketch.pnt.x[  iend] = xnew;
                sketch.pnt.y[  iend] = ynew;
                sketch.pnt.lbl[iend] = "";

                sketch.var.name[ nvar  ] = sprintf("::x[%]", npnt+1);
                sketch.var.value[nvar  ] = sprintf("%",    xnew  );

                sketch.var.name[ nvar+1] = sprintf("::y[%]", npnt+1);
                sketch.var.value[nvar+1] = sprintf("%",    ynew  );
            }

            // add the new Segment
            sketch.seg.type[nseg] = "B";
            sketch.seg.ibeg[nseg] = ibeg;
            sketch.seg.iend[nseg] = iend;
            sketch.seg.xm[  nseg] = (sketch.pnt.x[ibeg] + sketch.pnt.x[iend]) / 2;
            sketch.seg.ym[  nseg] = (sketch.pnt.y[ibeg] + sketch.pnt.y[iend]) / 2;
            sketch.seg.lbl[ nseg] = "";

            sketch.seg.dip[ nseg] = 0;
            sketch.seg.xc[  nseg] = 0;
            sketch.seg.yc[  nseg] = 0;
            sketch.seg.rad[ nseg] = 0;
            sketch.seg.tbeg[nseg] = 0;
            sketch.seg.tend[nseg] = 0;

            // add a new "horizontal" or "vertical" constraint
            if        (cnew == "H") {
                sketch.con.type[  ncon] = "H";
                sketch.con.index1[ncon] = ibeg;
                sketch.con.index2[ncon] = iend;
                sketch.con.value[ ncon] = 0;

                sketch.seg.lbl[nseg] += "H";
            } else if (cnew == "V") {
                sketch.con.type[  ncon] = "V";
                sketch.con.index1[ncon] = ibeg;
                sketch.con.index2[ncon] = iend;
                sketch.con.value[ ncon] = 0;

                sketch.seg.lbl[nseg] += "V";
            }

            // if Sketch is closed, post a message and swith mode
            if (iend == 0) {
                sketch.mode = 3;
            }

        // "z" --- add a zero-length Segment (to break a spline or Bezier)
        } else if (myKeyPress == 'z') {
            var ibeg = npnt - 1;
            var iend = 0;

            // make sure this follows a Bezier or spline
            if (sketch.seg.type[nseg-1] != "B" && sketch.seg.type[nseg-1] != "S") {
                postMessage("zero-length segment can only come after spline or bezier");
                return;
            }

            // add a new Point and its variables
            sketch.pnt.lbl[ibeg] = "Z";

            iend = npnt;
            sketch.pnt.x[  iend] = sketch.pnt.x[ibeg];
            sketch.pnt.y[  iend] = sketch.pnt.y[ibeg];
            sketch.pnt.lbl[iend] = "Z";

            sketch.var.name[ nvar  ] = sprintf("::x[%]", npnt);
            sketch.var.value[nvar  ] = sprintf("%",    xnew  );

            sketch.var.name[ nvar+1] = sprintf("::y[%]", npnt);
            sketch.var.value[nvar+1] = sprintf("%",    ynew  );

            // add a new (linear) Segment
            sketch.seg.type[nseg] = "L";
            sketch.seg.ibeg[nseg] = ibeg;
            sketch.seg.iend[nseg] = iend;
            sketch.seg.xm[  nseg] = sketch.pnt.x[ibeg];
            sketch.seg.ym[  nseg] = sketch.pnt.y[ibeg];
            sketch.seg.lbl[ nseg] = "Z";

            sketch.seg.dip[ nseg] = 0;
            sketch.seg.xc[  nseg] = 0;
            sketch.seg.yc[  nseg] = 0;
            sketch.seg.rad[ nseg] = 0;
            sketch.seg.tbeg[nseg] = 0;
            sketch.seg.tend[nseg] = 0;

            // add the constraints
            sketch.con.type[  ncon] = "Z";
            sketch.con.index1[ncon] = ibeg;
            sketch.con.index2[ncon] = -2;
            sketch.con.value[ ncon] = 0;

            sketch.con.type[  ncon+1] = "Z";
            sketch.con.index1[ncon+1] = iend;
            sketch.con.index2[ncon+1] = -3;
            sketch.con.value[ ncon+1] = 0;

        // "o" --- end (open) Sketch
        } else if (myKeyPress == "o") {
            sketch.mode = 3;

        } else {
            postMessage("Valid options are 'l', 'c', 's', 'b', 'z', and 'o'");
        }

    // sketch.mode 2 option (any keyPress or mouseDown)
    } else if (sketch.mode == 2) {

        // <any> --- set curvature
        if (sketch.seg.iend[sketch.seg.iend.length-1] == 0) {
            sketch.mode = 3;
        } else {
            sketch.mode = 1;
        }

    // sketch.mode 3 options
    // if at a Point:          'x', 'y', 'p', 't', 'a', 'w', 'd' or '<' to delete
    // if at a cirarc Segment: 'r', 's', 'x', 'y',               or '<' to delete
    // if at          Segment: 'h', 'v', 'i', 'l',               or '<' to delete
    } else if (sketch.mode == 3 || sketch.mode == 6) {
        sketcherSaveUndo();

        // "x" --- set X coordinate of nearest Point or arc Segment
        if (myKeyPress == "x" || myKeyPress == "X") {
            sketch.mode    = 3;
            sketch.suggest = undefined;

            var ibest = sketcherGetClosestPoint();
            var jbest = sketcherGetClosestSegment();
            if (jbest >= 0 && sketch.seg.type[jbest] != "C") {
                jbest = -1;
            }

            var dibest = 999999;
            if (ibest >= 0) {
                dibest = Math.abs(wv.cursorX - sketch.pnt.x[ibest])
                       + Math.abs(wv.cursorY - sketch.pnt.y[ibest]);
            }
            var djbest = 999999;
            if (jbest >= 0) {
                djbest = Math.abs(wv.cursorX - sketch.seg.xm[jbest])
                       + Math.abs(wv.cursorY - sketch.seg.ym[jbest]);
            }

            var xvalue;
            if (ibest < 0) {
                alert("No Point found");
            } else if (dibest < djbest) {
                if (sketch.pnt.lbl[ibest].indexOf("X") >= 0) {
                    for (var icon = 0; icon < ncon; icon++) {
                        if (sketch.con.type[icon] == "X" && sketch.con.index1[icon] == ibest) {
                            if (ibest == 0 && sketch.relative == 1) {
                                alert("X cannot be changed at first point if in relative mode");
                            } else {
                                xvalue = prompt("Enter x value: ", sketch.con.value[icon]);
                                if (xvalue !== null) {
                                    sketch.con.value[icon] = xvalue;
                                }
                            }
                            break;
                        }
                    }
                } else {
                    if (sketch.scale === undefined && sketch.xorig === undefined) {
                        xvalue = prompt("Enter x value: ");
                        if (xvalue !== null) {
                            for (var icon = 0; icon < ncon; icon++) {
                                if (sketch.con.type[icon] == "X" && sketch.con.index1[icon] != ibest) {
                                    sketch.scale = (xvalue - sketch.con.value[icon]) / (sketch.pnt.x[ibest] - sketch.pnt.x[sketch.con.index1[icon]]);
                                    sketch.xorig = sketch.pnt.x[sketch.con.index1[icon]] - sketch.con.value[icon] / sketch.scale;
                                    break;
                                }
                            }
                        }
                    } else if (sketch.xorig === undefined) {
                        xvalue = prompt("Enter x value: ");
                        if (xvalue !== null) {
                            sketch.xorig = sketch.pnt.x[ibest] + xvalue / sketch.scale;
                        }
                    } else if (sketch.scale === undefined) {
                        xvalue = prompt("Enter x value: ");
                        if (xvalue !== null) {
                            for (var icon = 0; icon < ncon; icon++) {
                                if (sketch.con.type[icon] == "X" && sketch.con.index1[icon] != ibest) {
                                    sketch.scale = (xvalue - sketch.con.value[icon]) / (sketch.pnt.x[ibest] - sketch.xorig);
                                    break;
                                }
                            }
                        }
                    } else {
                        xvalue = prompt("Enter x value: ", sketch.scale*(sketch.pnt.x[ibest]-sketch.xorig));
                    }

                    if (xvalue !== null) {
                        sketch.con.type[  ncon] = "X";
                        sketch.con.index1[ncon] = ibest;
                        sketch.con.index2[ncon] = -1;
                        sketch.con.value[ ncon] = xvalue;

                        sketch.pnt.lbl[ibest] += "X";
                    }
                }
            } else {
                if (sketch.seg.lbl[jbest].indexOf("X") >= 0) {
                    for (var icon = 0; icon < ncon; icon++) {
                        if (sketch.con.type[  icon] == "X" &&
                            sketch.con.index1[icon] == sketch.seg.ibeg[jbest] &&
                            sketch.con.index2[icon] == sketch.seg.iend[jbest]   ) {
                            xvalue = prompt("Enter x value: ", sketch.con.value[icon]);
                            if (xvalue !== null) {
                                sketch.con.value[icon] = xvalue;
                            }
                            break;
                        }
                    }
                } else {
                    xvalue = prompt("Enter x value: ");
                    if (xvalue !== null) {
                        sketch.con.type[  ncon] = "X";
                        sketch.con.index1[ncon] = sketch.seg.ibeg[jbest];
                        sketch.con.index2[ncon] = sketch.seg.iend[jbest];
                        sketch.con.value[ ncon] = xvalue;

                        sketch.seg.lbl[jbest] += "X";
                    }
                }
            }

        // "y" --- set Y coordinate of nearest Point or arc Segment
        } else if (myKeyPress == "y" || myKeyPress == "Y") {
            sketch.mode    = 3;
            sketch.suggest = undefined;

            var ibest = sketcherGetClosestPoint();
            var jbest = sketcherGetClosestSegment();
            if (jbest >= 0 && sketch.seg.type[jbest] != "C") {
                jbest = -1;
            }

            var dibest = 999999;
            if (ibest >= 0) {
                dibest = Math.abs(wv.cursorX - sketch.pnt.x[ibest])
                       + Math.abs(wv.cursorY - sketch.pnt.y[ibest]);
            }
            var djbest = 999999;
            if (jbest >= 0) {
                djbest = Math.abs(wv.cursorX - sketch.seg.xm[jbest])
                       + Math.abs(wv.cursorY - sketch.seg.ym[jbest]);
            }

            var yvalue;
            if (ibest < 0) {
                alert("No Point found");
            } else if (dibest < djbest) {
                if (sketch.pnt.lbl[ibest].indexOf("Y") >= 0) {
                    for (var icon = 0; icon < ncon; icon++) {
                        if (sketch.con.type[icon] == "Y" && sketch.con.index1[icon] == ibest) {
                            if (ibest == 0 && sketch.relative == 1) {
                                alert("Y cannot be changed at first point if in relative mode");
                            } else {
                                yvalue = prompt("Enter y value: ", sketch.con.value[icon]);
                                if (yvalue !== null) {
                                    sketch.con.value[icon] = yvalue;
                                }
                            }
                            break;
                        }
                    }
                } else {
                    if (sketch.scale === undefined && sketch.yorig === undefined) {
                        yvalue = prompt("Enter y value: ");
                        if (yvalue !== null) {
                            for (var icon = 0; icon < ncon; icon++) {
                                if (sketch.con.type[icon] == "Y" && sketch.con.index1[icon] != ibest) {
                                    sketch.scale = (sketch.con.value[icon] - yvalue) / (sketch.pnt.y[ibest] - sketch.pnt.y[sketch.con.index1[icon]]);
                                    sketch.yorig = sketch.pnt.y[sketch.con.index1[icon]] + sketch.con.value[icon] / sketch.scale;
                                    break;
                                }
                            }
                        }
                    } else if (sketch.yorig === undefined) {
                        yvalue = prompt("Enter y value: ");
                        if (yvalue !== null) {
                            sketch.yorig = sketch.pnt.y[ibest] - yvalue / sketch.scale;
                        }
                    } else if (sketch.scale === undefined) {
                        yvalue = prompt("Enter y value: ");
                        if (yvalue !== null) {
                            for (var icon = 0; icon < ncon; icon++) {
                                if (sketch.con.type[icon] == "Y" && sketch.con.index1[icon] != ibest) {
                                    sketch.scale = (sketch.con.value[icon] - yvalue) / (sketch.pnt.y[ibest] - sketch.yorig);
                                    break;
                                }
                            }
                        }
                    } else {
                        yvalue = prompt("Enter y value: ", sketch.scale*(sketch.yorig-sketch.pnt.y[ibest]));
                    }

                    if (yvalue !== null) {
                        sketch.con.type[  ncon] = "Y";
                        sketch.con.index1[ncon] = ibest;
                        sketch.con.index2[ncon] = -1;
                        sketch.con.value[ ncon] = yvalue;

                        sketch.pnt.lbl[ibest] += "Y";
                    }
                }
            } else {
                if (sketch.seg.lbl[jbest].indexOf("Y") >= 0) {
                    for (var icon = 0; icon < ncon; icon++) {
                        if (sketch.con.type[  icon] == "Y" &&
                            sketch.con.index1[icon] == sketch.seg.ibeg[jbest] &&
                            sketch.con.index2[icon] == sketch.seg.iend[jbest]   ) {
                            yvalue = prompt("Enter y value: ", sketch.con.value[icon]);
                            if (xvalue !== null) {
                                sketch.con.value[icon] = yvalue;
                            }
                            break;
                        }
                    }
                } else {
                    yvalue = prompt("Enter y value: ");
                    if (yvalue !== null) {
                        sketch.con.type[  ncon] = "Y";
                        sketch.con.index1[ncon] = sketch.seg.ibeg[jbest];
                        sketch.con.index2[ncon] = sketch.seg.iend[jbest];
                        sketch.con.value[ ncon] = yvalue;

                        sketch.seg.lbl[jbest] += "Y";
                    }
                }
            }

        // "p" --- perpendicularity
        } else if (myKeyPress == "p" || myKeyPress == "P") {
            sketch.mode    = 3;
            sketch.suggest = undefined;

            var ibest = sketcherGetClosestPoint();
            if (ibest < 0) {
                alert("No Point found");
            } else if (sketch.pnt.lbl[ibest].indexOf("P") >= 0) {
                alert("'P' constraint already at this Point");
            } else if (sketch.pnt.lbl[ibest].indexOf("T") >= 0) {
                alert("'T' constraint already at this Point");
            } else if (sketch.pnt.lbl[ibest].indexOf("A") >= 0) {
                alert("'A' constraint already at this Point");
            } else if (sketch.pnt.lbl[ibest].indexOf("Z") >= 0) {
                alert("'Z' constraint already at this Point");
            } else {
                sketch.con.type[  ncon] = "P";
                sketch.con.index1[ncon] = ibest;
                sketch.con.index2[ncon] = -1;
                sketch.con.value[ ncon] = 0;

                sketch.pnt.lbl[ibest] += "P";
            }

        // "t" --- tangency
        } else if (myKeyPress == "t" || myKeyPress == "T") {
            sketch.mode    = 3;
            sketch.suggest = undefined;

            var ibest = sketcherGetClosestPoint();
            if (ibest < 0) {
                alert("No Point found");
            } else if (sketch.pnt.lbl[ibest].indexOf("P") >= 0) {
                alert("'P' constraint already at this Point");
            } else if (sketch.pnt.lbl[ibest].indexOf("T") >= 0) {
                alert("'T' constraint already at this Point");
            } else if (sketch.pnt.lbl[ibest].indexOf("A") >= 0) {
                alert("'A' constraint already at this Point");
            } else if (sketch.pnt.lbl[ibest].indexOf("Z") >= 0) {
                alert("'Z' constraint already at this Point");
            } else {
                sketch.con.type[  ncon] = "T";
                sketch.con.index1[ncon] = ibest;
                sketch.con.index2[ncon] = -1;
                sketch.con.value[ ncon] = 0;

                sketch.pnt.lbl[ibest] += "T";
            }

        // "a" --- angle (positive to left)
        } else if (myKeyPress == 'a' || myKeyPress == "A") {
            sketch.mode    = 3;
            sketch.suggest = undefined;

            var ibest = sketcherGetClosestPoint();
            var avalue;
            if (ibest < 0) {
                alert("No Point found");
            } else if (sketch.pnt.lbl[ibest].indexOf("P") >= 0) {
                alert("'P' constraint already at this Point");
            } else if (sketch.pnt.lbl[ibest].indexOf("T") >= 0) {
                alert("'T' constraint already at this Point");
            } else if (sketch.pnt.lbl[ibest].indexOf("Z") >= 0) {
                alert("'Z' constraint already at this Point");
            } else if (sketch.pnt.lbl[ibest].indexOf("A") >= 0) {
                for (var icon = 0; icon < ncon; icon++) {
                    if (sketch.con.type[icon] == "A" && sketch.con.index1[icon] == ibest) {
                        avalue = prompt("Enter angle (deg): ", sketch.con.value[icon]);
                        if (avalue !== null) {
                            sketch.con.value[icon] = avalue;
                        }
                        break;
                    }
                }
            } else {
                var ibeg = -1;
                var iend = -1;
                for (var iseg = 0; iseg < sketch.seg.type.length; iseg++) {
                    if (sketch.seg.iend[iseg] == ibest) {
                        ibeg = sketch.seg.ibeg[iseg];
                    }
                    if (sketch.seg.ibeg[iseg] == ibest) {
                        iend = sketch.seg.iend[iseg];
                    }
                }

                if (ibeg < 0 || iend < 0) {
                    alert("Angle cannot be specified at endpoints");
                    return;
                }

                var ang1 = Math.atan2(sketch.pnt.y[ibest]-sketch.pnt.y[ibeg ],
                                      sketch.pnt.x[ibest]-sketch.pnt.x[ibeg ]) * 180 / Math.PI;
                var ang2 = Math.atan2(sketch.pnt.y[iend ]-sketch.pnt.y[ibest],
                                      sketch.pnt.x[iend ]-sketch.pnt.x[ibest]) * 180 / Math.PI;
                avalue = ang1 - ang2;
                while (avalue < -180) {
                    avalue += 360;
                }
                while (avalue > 180) {
                    avalue -= 360;
                }

                avalue = prompt("Enter angle (deg): ", avalue);

                if (avalue !== null) {
                    sketch.con.type[  ncon] = "A";
                    sketch.con.index1[ncon] = ibest;
                    sketch.con.index2[ncon] = -1;
                    sketch.con.value[ ncon] = avalue;

                    sketch.pnt.lbl[ibest] += "A";
                }
            }

        // "h" --- Segment is horizontal
        } else if (myKeyPress == "h" || myKeyPress == "H") {
            sketch.mode    = 3;
            sketch.suggest = undefined;

            var ibest = sketcherGetClosestSegment();
            if (ibest < 0) {
                alert("No Segment found");
            } else if (sketch.seg.lbl[ibest].indexOf("H") >= 0) {
                alert("'H' constraint already on this Segment");
            } else if (sketch.seg.lbl[ibest].indexOf("V") >= 0) {
                alert("'V' constraint already on this Segment");
            } else if (sketch.seg.lbl[ibest].indexOf("I") >= 0) {
                alert("'I' constraint already on this Segment");
            } else if (Math.abs(sketch.seg.dip[ibest]) > 1e-3) {
                alert("'H' constraint cannot be set on a Cirarc");
            } else {
                var ibeg  = sketch.seg.ibeg[ibest];
                var iend  = sketch.seg.iend[ibest];

                sketch.con.type[  ncon] = "H";
                sketch.con.index1[ncon] = ibeg;
                sketch.con.index2[ncon] = iend;
                sketch.con.value[ ncon] = 0;

                sketch.seg.lbl[ibest] += "H";
            }

        // "v" --- Segment is vertical
        } else if (myKeyPress == "v" || myKeyPress == "V") {
            sketch.mode    = 3;
            sketch.suggest = undefined;

            var ibest = sketcherGetClosestSegment();
            if (ibest < 0) {
                alert("No Segment found");
            } else if (sketch.seg.lbl[ibest].indexOf("V") >= 0) {
                alert("'V' constraint already on this Segment");
            } else if (sketch.seg.lbl[ibest].indexOf("H") >= 0) {
                alert("'H' constraint already on this Segment");
            } else if (sketch.seg.lbl[ibest].indexOf("I") >= 0) {
                alert("'I' constraint already on this Segment");
            } else if (Math.abs(sketch.seg.dip[ibest]) > 1e-3) {
                alert("'V' constraint cannot be set on a Cirarc");
            } else {
                var ibeg  = sketch.seg.ibeg[ibest];
                var iend  = sketch.seg.iend[ibest];

                sketch.con.type[  ncon] = "V";
                sketch.con.index1[ncon] = ibeg;
                sketch.con.index2[ncon] = iend;
                sketch.con.value[ ncon] = 0;

                sketch.seg.lbl[ibest] += "V";
            }

        // "i" --- inclination of Segment
        } else if (myKeyPress == "i" || myKeyPress == "I") {
            sketch.mode    = 3;
            sketch.suggest = undefined;

            var ibest = sketcherGetClosestSegment();
            var ivalue;
            if (ibest < 0) {
                alert("No Segment found");
            } else if (sketch.seg.lbl[ibest].indexOf("H") >= 0) {
                alert("'H' constraint already on this Segment");
            } else if (sketch.seg.lbl[ibest].indexOf("V") >= 0) {
                alert("'V' constraint already on this Segment");
            } else if (Math.abs(sketch.seg.dip[ibest]) > 1e-3) {
                alert("'I' constraint cannot be set on a Cirarc");
            } else if (sketch.seg.lbl[ibest].indexOf("I") >= 0) {
                for (var icon = 0; icon < ncon; icon++) {
                    if (sketch.con.type[ icon] == "I" && sketch.con.index1[icon] == ibest) {
                        ivalue = prompt("Enter inclination angle (deg): ", sketch.con.value[icon]);
                        if (ivalue !== null) {
                            sketch.con.value[icon] = ivalue;
                        }
                        break;
                    }
                }
            } else {
                var ibeg  = sketch.seg.ibeg[ibest];
                var iend  = sketch.seg.iend[ibest];

                var iprompt = Math.atan2(sketch.pnt.y[ibeg]-sketch.pnt.y[iend], sketch.pnt.x[iend]-sketch.pnt.x[ibeg]) * 180/Math.PI;

                ivalue = prompt("Enter inclination angle (deg): ", iprompt);
                if (ivalue !== null) {
                    sketch.con.type[  ncon] = "I";
                    sketch.con.index1[ncon] = ibeg;
                    sketch.con.index2[ncon] = iend;
                    sketch.con.value[ ncon] = ivalue;

                    sketch.seg.lbl[ibest] += "I";
                }
            }

        // "l" --- length of Segment
        } else if (myKeyPress == "l" || myKeyPress == "L") {
            sketch.mode    = 3;
            sketch.suggest = undefined;

            var ibest = sketcherGetClosestSegment();
            var lvalue;
            if (ibest < 0) {
                alert("No Segment found");
            } else if (Math.abs(sketch.seg.dip[ibest]) > 1e-3) {
                alert("'L' constraint cannot be set on a Cirarc");
            } else if (sketch.seg.lbl[ibest].indexOf("L") >= 0) {
                for (var icon = 0; icon < ncon; icon++) {
                    if (sketch.con.type[ icon] == "L" && sketch.con.index1[icon] == ibest) {
                        lvalue = prompt("Enter length: ", sketch.con.value[icon]);
                        if (lvalue !== null) {
                            sketch.con.value[icon] = lvalue;
                        }
                        break;
                    }
                }
            } else {
                var ibeg  = sketch.seg.ibeg[ibest];
                var iend  = sketch.seg.iend[ibest];

                var xbeg = sketch.pnt.x[ibeg];
                var ybeg = sketch.pnt.y[ibeg];
                var xend = sketch.pnt.x[iend];
                var yend = sketch.pnt.y[iend];
                var len  = Math.sqrt(Math.pow(xend-xbeg, 2) + Math.pow(yend-ybeg, 2));

                if (sketch.scale === undefined) {
                    lvalue = prompt("Enter length: ");
                    if (lvalue !== null) {
                        sketch.scale = lvalue / len;
                    }
                } else {
                    lvalue = prompt("Enter length: ", sketch.scale*len);
                }

                if (lvalue !== null) {
                    sketch.con.type[  ncon] = "L";
                    sketch.con.index1[ncon] = ibeg;
                    sketch.con.index2[ncon] = iend;
                    sketch.con.value[ ncon] = lvalue;

                    sketch.seg.lbl[ibest] += "L";
                }
            }

        // "r" --- radius of Segment
        } else if (myKeyPress == "r" || myKeyPress == "R") {
            sketch.mode    = 3;
            sketch.suggest = undefined;

            var ibest = sketcherGetClosestSegment();
            var rvalue;
            if (ibest < 0) {
                alert("No Cirarc found");
            } else if (Math.abs(sketch.seg.dip[ibest]) < 1e-3) {
                alert("'R' constraint cannot be set for a Linseg");
            } else if (sketch.seg.lbl[ibest].indexOf("R") >= 0) {
                for (var icon = 0; icon < ncon; icon++) {
                    if (sketch.con.type[ icon] == "R" && sketch.con.index1[icon] == ibest) {
                        rvalue = prompt("Enter radius: ", sketch.con.value[icon]);
                        if (rvalue !== null) {
                            sketch.con.value[icon] = rvalue;
                        }
                        break;
                    }
                }
            } else {
                var ibeg = sketch.seg.ibeg[ibest];
                var iend = sketch.seg.iend[ibest];

                var xa = sketch.pnt.x[ibeg];
                var ya = sketch.pnt.y[ibeg];
                var xb = sketch.pnt.x[iend];
                var yb = sketch.pnt.y[iend];

                var d  = sketch.seg.dip[ibest];

                var L   = Math.sqrt((xb-xa) * (xb-xa) + (yb-ya) * (yb-ya));
                var rad = (L*L + 4*d*d) / (8*d);

                if (sketch.scale === undefined) {
                    if (Number(sketch.seg.dip[ibest]) >= 0) {
                        rvalue = prompt("Enter radius (should be positive as drawn): ");
                    } else {
                        rvalue = prompt("Enter radius (should be negative as drawn): ");
                    }
                    if (rvalue !== null) {
                        sketch.scale = rvalue / rad;
                    }
                } else {
                    rvalue = prompt("Enter radius: ", sketch.scale*rad);
                }

                if (rvalue !== null) {
                    sketch.con.type[  ncon] = "R";
                    sketch.con.index1[ncon] = ibeg;
                    sketch.con.index2[ncon] = iend;
                    sketch.con.value[ ncon] = rvalue;

                    sketch.seg.lbl[ibest] += "R";
                }
            }

        // "s" --- sweep angle of Segment
        } else if (myKeyPress == "s" || myKeyPress == "S") {
            sketch.mode    = 3;
            sketch.suggest = undefined;

            var ibest = sketcherGetClosestSegment();
            var avalue;
            if (ibest < 0) {
                alert("No Cirarc found");
            } else if (Math.abs(sketch.seg.dip[ibest]) < 1e-3) {
                alert("'S' constraint cannot be set for a Linseg");
            } else if (sketch.seg.lbl[ibest].indexOf("S") >= 0) {
                for (var icon = 0; icon < ncon; icon++) {
                    if (sketch.con.type[icon] == "S" && sketch.con.index1[icon] == ibest) {
                        avalue = prompt("Enter sweep angle (deg): ", sketch.con.value[icon]);
                        if (avalue !== null) {
                            sketch.con.value[icon] = avalue;
                        }
                        break;
                    }
                }
            } else {
                var ibeg = sketch.seg.ibeg[ibest];
                var iend = sketch.seg.iend[ibest];

                if (Number(sketch.seg.dip[ibest]) >= 0) {
                    avalue = prompt("Enter sweep angle (deg, should be positive as drawn): ");
                } else {
                    avalue = prompt("Enter sweep angle (deg, should be negative as drawn): ");
                }
                if (avalue !== null) {
                    sketch.con.type[  ncon] = "S";
                    sketch.con.index1[ncon] = ibeg;
                    sketch.con.index2[ncon] = iend;
                    sketch.con.value[ ncon] = avalue;

                    sketch.seg.lbl[ibest] += "S";
                }
            }

        // "w" --- width (dx) between Points
        } else if (myKeyPress == "w" || myKeyPress == "W") {
            var bestPnt = sketcherGetClosestPoint();
            var bestCon = sketcherGetClosestConstraint();

            if (bestCon !== undefined) {
                var xmid = (  sketch.pnt.x[sketch.con.index1[bestCon]]
                            + sketch.pnt.x[sketch.con.index2[bestCon]]) / 2;
                var ymid = (2*sketch.pnt.y[sketch.con.index1[bestCon]]
                            + sketch.pnt.y[sketch.con.index2[bestCon]]) / 3;

                var dpnt = Math.abs(sketch.pnt.x[bestPnt] - wv.cursorX)
                         + Math.abs(sketch.pnt.y[bestPnt] - wv.cursorY);
                var dcon = Math.abs(xmid                - wv.cursorX)
                         + Math.abs(ymid                - wv.cursorY);

                if (dcon < dpnt) {
                    var wvalue = prompt("Enter width: ", sketch.con.value[bestCon]);
                    if (wvalue !== null) {
                        sketch.con.value[bestCon] = wvalue;
                    }

                    sketch.bestPoint = undefined;
                    sketch.mode      = 3;
                    sketch.suggest   = undefined;
                } else {
                    sketch.basePoint = bestPnt;
                    sketch.mode      = 4;
                    sketch.suggest   = undefined;
                }
            } else {
                sketch.basePoint = bestPnt;
                sketch.mode      = 4;
                sketch.suggest   = undefined;
            }

        // "d" --- depth (dy) between Points
        } else if (myKeyPress == "d" || myKeyPress == "D") {
            var bestPnt = sketcherGetClosestPoint();
            var bestCon = sketcherGetClosestConstraint();

            if (bestCon !== undefined) {
                var xmid = (2*sketch.pnt.x[sketch.con.index1[bestCon]]
                            + sketch.pnt.x[sketch.con.index2[bestCon]]) / 3;
                var ymid = (  sketch.pnt.y[sketch.con.index1[bestCon]]
                            + sketch.pnt.y[sketch.con.index2[bestCon]]) / 2;

                var dpnt = Math.abs(sketch.pnt.x[bestPnt] - wv.cursorX)
                         + Math.abs(sketch.pnt.y[bestPnt] - wv.cursorY);
                var dcon = Math.abs(xmid                - wv.cursorX)
                         + Math.abs(ymid                - wv.cursorY);

                if (dcon < dpnt) {
                    var dvalue = prompt("Enter depth: ", sketch.con.value[bestCon]);
                    if (dvalue !== null) {
                        sketch.con.value[bestCon] = dvalue;
                    }

                    sketch.bestPoint = undefined;
                    sketch.mode      = 3;
                    sketch.suggest   = undefined;
                } else {
                    sketch.basePoint = bestPnt;
                    sketch.mode      = 5;
                    sketch.suggest   = undefined;
                }
            } else {
                sketch.basePoint = bestPnt;
                sketch.mode      = 5;
                sketch.suggest   = undefined;
            }

        // "<" --- remove constraint(s) from nearset Point or Segment
        } else if (myKeyPress == "<") {
            sketch.mode    = 3;
            sketch.suggest = undefined;

            var ibest = sketcherGetClosestPoint();
            var jbest = sketcherGetClosestSegment();
            var kbest = sketcherGetClosestConstraint();

            var dibest = 999999;
            if (ibest >= 0) {
                dibest = Math.abs(wv.cursorX - sketch.pnt.x[ibest])
                       + Math.abs(wv.cursorY - sketch.pnt.y[ibest]);
            }
            var djbest = 999999;
            if (jbest >= 0) {
                djbest = Math.abs(wv.cursorX - sketch.seg.xm[jbest])
                       + Math.abs(wv.cursorY - sketch.seg.ym[jbest]);
            }
            var dkbest = 999999;
            if (kbest !== undefined) {
                if        (sketch.con.type[kbest] == "W") {
                    var xmid = (  sketch.pnt.x[sketch.con.index1[kbest]]
                                + sketch.pnt.x[sketch.con.index2[kbest]]) / 2;
                    var ymid = (2*sketch.pnt.y[sketch.con.index1[kbest]]
                                + sketch.pnt.y[sketch.con.index2[kbest]]) / 3;

                    dkbest = Math.abs(wv.cursorX - xmid)
                           + Math.abs(wv.cursorY - ymid);
                } else if (sketch.con.type[kbest] == "D") {
                    var xmid = (2*sketch.pnt.x[sketch.con.index1[kbest]]
                                + sketch.pnt.x[sketch.con.index2[kbest]]) / 3;
                    var ymid = (  sketch.pnt.y[sketch.con.index1[kbest]]
                                + sketch.pnt.y[sketch.con.index2[kbest]]) / 2;

                    dkbest = Math.abs(wv.cursorX - xmid)
                           + Math.abs(wv.cursorY - ymid);
                }
            }

            // remove (selected) constraint at Point ibest
            if (ibest >= 0 && dibest < djbest && dibest < dkbest) {
                if (sketch.pnt.lbl[ibest].length == 0) {
                    postMessage("no constraint at Point "+ibest);
                } else if (sketch.pnt.lbl[ibest] == "Z") {
                    postMessage("cannot delete a 'Z' constraint");
                } else {
                    var jtype;
                    if (sketch.pnt.lbl[ibest].length == 1) {
                        jtype = sketch.pnt.lbl[ibest];
                    } else {
                        jtype = prompt("More than one constraint at Point "+ibest+
                                      "\nEnter one character from \""+sketch.pnt.lbl[ibest]+
                                      "\" to remove");
                        if (jtype === null) {
                            jtype = "&";
                        }
                    }
                    jtype = jtype.toUpperCase();

                    if        (ibest == 0 && sketch.relative == 1 && jtype == "X") {
                        alert("X cannot be deleted from first point if in relative mode");
                    } else if (ibest == 0 && sketch.relative == 1 && jtype == "Y") {
                        alert("Y cannot be deleted from first point if in relative mode");
                    } else {
                        for (var icon = 0; icon < sketch.con.type.length; icon++) {
                            if (sketch.con.index1[icon] == ibest &&
                                sketch.con.index2[icon] <  0     &&
                                sketch.con.type[  icon] == jtype   ) {
                                postMessage("removing \""+jtype+"\" contraint from Point "+ibest);

                                sketch.con.type.splice(  icon, 1);
                                sketch.con.index1.splice(icon, 1);
                                sketch.con.index2.splice(icon, 1);
                                sketch.con.value.splice( icon, 1);

                                sketch.pnt.lbl[ibest] = sketch.pnt.lbl[ibest].split(jtype).join("");
                                break;
                            }
                        }
                    }
                }

            // remove constraints at Segment jbest
            } else if (jbest >= 0 && djbest < dkbest) {
                if (sketch.seg.lbl[jbest].length == 0) {
                    postMessage("no constraint at Segment "+jbest);
                } else if (sketch.seg.lbl[ibest] == "Z") {
                    postMessage("cannot delete a 'Z' constraint");
                } else {
                    var jtype;
                    if (sketch.seg.lbl[jbest].length == 1) {
                        jtype = sketch.seg.lbl[jbest];
                    } else {
                        jtype = prompt("More than one constraint at Segment "+jbest+
                                      "\nEnter one character from \""+sketch.seg.lbl[jbest]+
                                      "\" to remove");
                        if (jtype === null) {
                            jtype = "&";
                        }
                    }
                    jtype = jtype.toUpperCase();

                    for (var icon = 0; icon < sketch.con.type.length; icon++) {
                        if (sketch.con.index1[icon] == jbest &&
                            sketch.con.index2[icon] >= 0     &&
                            sketch.con.type[  icon] == jtype   ) {
                            postMessage("removing \""+jtype+"\" constraint from Segment "+jbest);

                            sketch.con.type.splice(  icon, 1);
                            sketch.con.index1.splice(icon, 1);
                            sketch.con.index2.splice(icon, 1);
                            sketch.con.value.splice( icon, 1);

                            sketch.seg.lbl[jbest] = sketch.seg.lbl[jbest].split(jtype).join("");
                            break;
                        }
                    }
                }

            // remove H or W constraint
            } else {
                postMessage("removing constraint "+kbest);

                sketch.con.type.splice(  kbest, 1);
                sketch.con.index1.splice(kbest, 1);
                sketch.con.index2.splice(kbest, 1);
                sketch.con.value.splice( kbest, 1);
            }

        // "?" --- examine Point or Segment
        } else if (myKeyPress == "?") {
            var ibest = sketcherGetClosestPoint();
            var jbest = sketcherGetClosestSegment();

            var dibest = 999999;
            if (ibest >= 0) {
                dibest = Math.abs(wv.cursorX - sketch.pnt.x[ibest])
                       + Math.abs(wv.cursorY - sketch.pnt.y[ibest]);
            }
            var djbest = 999999;
            if (jbest >= 0) {
                djbest = Math.abs(wv.cursorX - sketch.seg.xm[jbest])
                       + Math.abs(wv.cursorY - sketch.seg.ym[jbest]);
            }

            // examine Point ibest
            if (ibest >= 0 && dibest < djbest) {
                if (sketch.pnt.lbl[ibest].length == 0) {
                    postMessage("Point "+(ibest+1)+" has no constraints");
                } else {
                    postMessage("Point "+(ibest+1)+" has constraints: "+sketch.pnt.lbl[ibest]);
                }

            // examine Segment jbest
            } else if (jbest >= 0) {
                if        (sketch.seg.type[jbest] == "L") {
                    var segtype = "linseg";
                } else if (sketch.seg.type[jbest] == "C") {
                    var segtype = "cirarc";
                } else if (sketch.seg.type[jbest] == "S") {
                    var segtype = "spline";
                } else if (sketch.seg.type[jbest] == "B") {
                    var segtype = "bezier";
                } else {
                    var segtype = "**unknown type**";
                }
                if (sketch.seg.lbl[jbest].length == 0) {
                    postMessage("Segment "+(jbest+1)+" ("+segtype+") has no constraints");
                } else {
                    postMessage("Segment "+(jbest+1)+" ("+segtype+") has constraints: "+sketch.seg.lbl[jbest]);
                }
            }

        // "@" --- dump sket structure
        } else if (myKeyPress == "@") {
            postMessage("sketch.pnt->"+JSON.stringify(sketch.pnt));
            postMessage("sketch.seg->"+JSON.stringify(sketch.seg));
            postMessage("sketch.var->"+JSON.stringify(sketch.var));
            postMessage("sketch.con->"+JSON.stringify(sketch.con));

        // invalid option
        } else {
            var ibest = sketcherGetClosestPoint();
            var jbest = sketcherGetClosestSegment();

            var dibest = 999999;
            if (ibest >= 0) {
                dibest = Math.abs(wv.cursorX - sketch.pnt.x[ibest])
                       + Math.abs(wv.cursorY - sketch.pnt.y[ibest]);
            }
            var djbest = 999999;
            if (jbest >= 0) {
                djbest = Math.abs(wv.cursorX - sketch.seg.xm[jbest])
                       + Math.abs(wv.cursorY - sketch.seg.ym[jbest]);
            }

            // remove constraints at Point ibest
            if (ibest >= 0 && dibest < djbest) {
                postMessage("Valid options are 'x', 'y', 'p', 't', 'w', 'd', and '<'");
            } else {
                postMessage("Valid options are 'h', 'v', 'i', 'l', 'r', 's', and '<'");
            }
        }

        // if .scale has been set, try to update .xorig or .yorig
        //    if a suitable constraint exists
        if (sketch.scale !== undefined) {
            if (sketch.xorig === undefined) {
                for (var icon = 0; icon < ncon; icon++) {
                    if (sketch.con.type[icon] == "X") {
                        sketch.xorig = sketch.pnt.x[sketch.con.index1[icon]] - sketch.con.value[icon] / sketch.scale;
                        break;
                    }
                }
            }
            if (sketch.yorig === undefined) {
                for (var icon = 0; icon < ncon; icon++) {
                    if (sketch.con.type[icon] == "Y") {
                        sketch.yorig = sketch.pnt.y[sketch.con.index1[icon]] + sketch.con.value[icon] / sketch.scale;
                        break;
                    }
                }
            }
        }

    // sketch.mode 4 options: w
    } else if (sketch.mode == 4) {
        sketch.mode    = 3;
        sketch.suggest = undefined;

        var target = sketcherGetClosestPoint();

        if (target == sketch.basePoint) {
            alert("Width not set since base and target Point are the same");
        } else {
            var wvalue = (sketch.pnt.x[target] - sketch.pnt.y[sketch.basePoint]) * sketch.scale;
            if (wvalue > 0) {
                wvalue = prompt("Enter width (should be positive as drawn): ", wvalue);
            } else {
                wvalue = prompt("Enter width (should be negative as drawn): ", wvalue);
            }
            if (wvalue !== null) {
                sketch.con.type[  ncon] = "W";
                sketch.con.index1[ncon] = sketch.basePoint;
                sketch.con.index2[ncon] = target;
                sketch.con.value[ ncon] = wvalue;
            }
        }

        sketch.basePoint = undefined;

    // sketch.mode 5 options: d
    } else if (sketch.mode == 5) {
        sketch.mode    = 3;
        sketch.suggest = undefined;

        var target = sketcherGetClosestPoint();

        if (target == sketch.basePoint) {
            alert("Depth not set since base and target Point are the same");
        } else {
            var dvalue = (sketch.pnt.y[sketch.basePoint] - sketch.pnt.y[target]) * sketch.scale;
            if (dvalue > 0) {
                dvalue = prompt("Enter depth (should be positive as drawn): ", dvalue);
            } else {
                dvalue = prompt("Enter depth (should be negative as drawn): ", dvalue);
            }
            if (dvalue !== null) {
                sketch.con.type[  ncon] = "D";
                sketch.con.index1[ncon] = sketch.basePoint;
                sketch.con.index2[ncon] = target;
                sketch.con.value[ ncon] = dvalue;
            }
        }

        sketch.basePoint = undefined;
    }

    sketcherDraw();
};

////////////////////////////////////////////////////////////////////////////////


//
// Initialize a Sketch
//
var sketcherInitialize = function () {
    // alert("sketcherInitialize()");

    sketch.basePoint = undefined;
    sketch.skbegX    = undefined;
    sketch.skbegY    = undefined;
    sketch.skbegZ    = undefined;
    sketch.relative  =  0;
    sketch.suggest   = undefined;
    sketch.heads     = [];
    sketch.tails     = [];
    sketch.halo      =  5;      // pixels for determining "closeness" while drawing
    wv.keyPress    = -1;      // last key pressed
    wv.keyCode     = -1;
    wv.button      = -1;      // button pressed
    wv.modifier    =  0;      // modifier (shift,alt,ctrl) bitflag

    // set up Sketcher offsets
    var sketcher = document.getElementById("sketcher");
    sketch.offTop    = sketcher.offsetTop;
    sketch.offLeft   = sketcher.offsetLeft;

    // initialize the arrays
    sketch.pnt        = {};
    sketch.pnt.x      = [];
    sketch.pnt.y      = [];
    sketch.pnt.lbl    = [];

    sketch.seg        = {};
    sketch.seg.type   = [];
    sketch.seg.ibeg   = [];
    sketch.seg.iend   = [];
    sketch.seg.xm     = [];
    sketch.seg.ym     = [];
    sketch.seg.lbl    = [];

    sketch.seg.dip    = [];
    sketch.seg.xc     = [];
    sketch.seg.yc     = [];
    sketch.seg.rad    = [];
    sketch.seg.tbeg   = [];
    sketch.seg.tend   = [];

    sketch.var        = {};
    sketch.var.name   = [];
    sketch.var.value  = [];

    sketch.con        = {};
    sketch.con.type   = [];
    sketch.con.index1 = [];
    sketch.con.index2 = [];
    sketch.con.value  = [];

    // mapping between "physical" and "screen" coordinates
    //   x_physical = sketch.scale * (x_screen - sketch.xorig)
    //   y_physical = sketch.scale * (sketch.yorig - y_screen)

    // commands that set sketch.xorig: "x"
    // commands that set sketch.yorig: "y"
    // commands that set sketch.scale: "l", "r", "w", or "d"
    //                               "x" if another "x" exists
    //                               "y" if another "y" exists

    // reset scaling info
    sketch.scale = undefined;
    sketch.xorig = undefined;
    sketch.yorig = undefined;

    // initialize undo
    sketch.maxundo = 10;
    sketch.nundo   =  0;
    sketch.undo    = new Array;
    for (var iundo = 0; iundo < sketch.maxundo; iundo++) {
        sketch.undo[iundo] = undefined;
    }

    // set mode to initializing
    sketch.mode = 0;

    // change done button legend
    var button = document.getElementById("doneMenuBtn");
    button["innerHTML"] = "Sketch";

    // draw the current Sketch
    sketcherDraw();
};


//
// called when "solveButton" is pressed
//
var sketcherSolvePre = function () {
    // alert("sketcherSolvePre()");

    var npnt = sketch.pnt.x.length;
    var nvar = sketch.var.name.length;
    var ncon = sketch.con.type.length;
    var nseg = sketch.seg.type.length;

    // if (nvar > ncon) {
    //     alert("Sketch is under-constrained since nvar="+nvar+" but ncon="+ncon);
    //     return;
    // } else if (nvar < ncon) {
    //     alert("Sketch is over-constrained since nvar="+nvar+" but ncon="+ncon);
    //     return;
    // }

    if (sketch.xorig === undefined || isNaN(sketch.xorig) ||
        sketch.yorig === undefined || isNaN(sketch.yorig) ||
        sketch.scale === undefined || isNaN(sketch.scale)   ) {

        // get the size of the canvas
        var canvas = document.getElementById("sketcher");
        var width  = canvas.clientWidth;
        var height = canvas.clientHeight;

        // guess that scale of body is "5" and that it is centered at origin
        sketch.scale = 5 / width;
        sketch.xorig = width  / 2;
        sketch.yorig = height / 2;
    }

    // create a message to send to the server
    var message = "solveSketch|";

    for (var ipnt = 0; ipnt < npnt; ipnt++) {
        if (ipnt == 0) {
            var im1 = npnt - 1;
        } else {
            var im1 = ipnt - 1;
        }

        message += Number(sketch.scale * (sketch.pnt.x[ipnt] - sketch.xorig)).toFixed(6) + ";";
        message += Number(sketch.scale * (sketch.yorig - sketch.pnt.y[ipnt])).toFixed(6) + ";";

        if (ipnt > 0) {
            message += Number(sketch.scale *  sketch.seg.dip[ipnt-1]      ).toFixed(6) + ";";
        } else if (sketch.seg.iend[nseg-1] == 0) {
            message += Number(sketch.scale *  sketch.seg.dip[npnt-1]      ).toFixed(6) + ";";
        } else {
            message += "0;";
        }
    }

    message += "|";

    for (var icon = 0; icon < ncon; icon++) {
        message += sketch.con.type[icon] + ";" + (sketch.con.index1[icon]+1) + ";";
        if (sketch.con.index2[icon] >= 0) {
            message += (sketch.con.index2[icon] + 1) + ";" + sketch.con.value[icon] + ";";
        } else {
            message +=  sketch.con.index2[icon]      + ";" + sketch.con.value[icon] + ";";
        }
    }

    message += "|";

    browserToServer(message);

//    // change button legend
//    var button = document.getElementById("solveButton");
//    button["innerHTML"] = "Solving...";
//    button.style.backgroundColor = "#FFFF3F";
};


//
// update the variables associated with a Sketch
//
var sketcherSolvePost = function (text) {
    // alert("sketcherSolvePost(text="+text+")");

    // break up based upon the "|"
    var textList = text.split("|");

    if        (textList[1].substring(0,7) == "ERROR::") {
        postMessage(textList[1]);
        alert("Sketch was not successfully solved, but no hints are available\n" +
              "Note: try adding a R or S constraint to a circular arc.");
        return;
    } else if (textList[1].substring(0,4) == "*del") {
        sketch.suggest = textList[1];
        sketcherDraw();
        postMessage("Delete one of the constraints in red (using < key)\n" +
                    "Note: others may be possible as well.");
        return;
    } else if (textList[1].substring(0,4) == "*add") {
        sketch.suggest = textList[1];
        sketcherDraw();
        postMessage("Add one of the constraints in green\n" +
                    "Note: others may be possible as well.");
        return;
    }

    var npnt = 0;
    for (var i = 0; i < textList[1].length; i++) {
        if (textList[1].charAt(i) == ';') {
            npnt++;
        }
    }
    npnt /= 3;

    // update the values (which are separated by semicolons)
    var values = textList[1].split(";");

    for (var ipnt = 0; ipnt < npnt; ipnt++) {
        if (ipnt == 0) {
            var im1 = npnt - 1;
        } else {
            var im1 = ipnt - 1;
        }

//        sketch.pnt.x[  ipnt] = Math.floor(sketch.xorig + Number(values[3*ipnt  ]) / sketch.scale);
//        sketch.pnt.y[  ipnt] = Math.floor(sketch.yorig - Number(values[3*ipnt+1]) / sketch.scale);
//        sketch.seg.dip[im1 ] = Math.floor(             Number(values[3*ipnt+2]) / sketch.scale);
        sketch.pnt.x[  ipnt] = sketch.xorig + Number(values[3*ipnt  ]) / sketch.scale;
        sketch.pnt.y[  ipnt] = sketch.yorig - Number(values[3*ipnt+1]) / sketch.scale;
        sketch.seg.dip[im1 ] =              Number(values[3*ipnt+2]) / sketch.scale;
    }

    // update data associated with segments
    for (var iseg = 0; iseg < sketch.seg.type.length; iseg++) {
        sketcherUpdateSegData(iseg);
    }

    // set the mode
    sketch.mode = 6;

    // redraw the Sketch with the updated values
    sketcherDraw();
};


//
// called when "save" is pressed
//
var sketcherSave = function () {
    // alert("in sketcherSave()");

    // if the scaling is not set, make sure that at least one length is given
    //    (or else NaN will be sent back to .csm)
    if        (sketch.scale === undefined || isNaN(sketch.scale)) {
        alert("At lest one length or radius must be set before saving");
        return;
    } else if (sketch.xorig === undefined || isNaN(sketch.xorig)) {
        alert("At least one \"X\" must be set before saving");
        return;
    } else if (sketch.yorig === undefined || isNaN(sketch.yorig)) {
        alert("At least one \"Y\" must be set before saving");
        return;
    }

    // alert user that an unsolved Sketch may cause .csm to not solve
    if (sketch.mode != 6) {
        if (confirm("Sketch is not solved and may break build.  Continue?") !== true) {
            return;
        }
    }

    var npnt = sketch.pnt.x.length;
    var nseg = sketch.seg.type.length;
    var ncon = sketch.con.type.length;

    // build a string of the variables
    var vars = "";
    for (var ipnt = 0; ipnt < npnt; ipnt++) {
        vars += Number(sketch.scale * (sketch.pnt.x[ipnt] - sketch.xorig)).toFixed(6) + ";";
        vars += Number(sketch.scale * (sketch.yorig - sketch.pnt.y[ipnt])).toFixed(6) + ";";

        if (ipnt > 0) {
            vars += Number(sketch.scale *  sketch.seg.dip[ipnt-1]      ).toFixed(6) + ";";
        } else if (sketch.seg.iend[nseg-1] == 0) {
            vars += Number(sketch.scale *  sketch.seg.dip[npnt-1]      ).toFixed(6) + ";";
        } else {
            vars += "0;";
        }
    }

    // build a string of the constraints
    var cons = "";
    for (var icon = 0; icon < ncon; icon++) {
        cons += sketch.con.type[icon] + ";" + (sketch.con.index1[icon]+1) + ";";
        if (sketch.con.index2[icon] >= 0) {
            cons += (sketch.con.index2[icon] + 1) + ";" + sketch.con.value[icon] + ";";
        } else {
            cons +=  sketch.con.index2[icon]      + ";" + sketch.con.value[icon] + ";";
        }
    }
    cons = cons.replace(/\s/g, "");

    // build a string of the segments
    var segs = "";
    for (var iseg = 0; iseg < nseg; iseg++) {
        segs += sketch.seg.type[iseg] + ";" + (sketch.seg.ibeg[iseg]+1) + ";" +
                                            (sketch.seg.iend[iseg]+1) + ";";
    }

    // create the message to be sent to the server
    var mesg = ""+sketch.ibrch+"|"+vars+"|"+cons+"|"+segs+"|";

    // because of an apparent limit on the size of text
    //    messages that can be sent from the browser to the
    //    server, we need to send the new file back in
    //    pieces and then reassemble on the server
    var maxMessageSize = 1000;

    var ichar = 0;
    var part  = mesg.substring(ichar, ichar+maxMessageSize);
    browserToServer("saveSketchBeg|"+part);
    ichar += maxMessageSize;

    while (ichar < mesg.length) {
        part = mesg.substring(ichar, ichar+maxMessageSize);
        browserToServer("saveSketchMid|"+part);
        ichar += maxMessageSize;
    }

    browserToServer("saveSketchEnd|");
};


//
// undo last Sketch command
//
var sketcherUndo = function () {
    // alert("sketcherUndo()");

    var iundo = (sketch.nundo-1) % sketch.maxundo;

    // if undo information does not exist, inform the user
    if (sketch.undo[iundo]     === undefined ||
        sketch.undo[iundo].pnt === undefined ||
        sketch.undo[iundo].seg === undefined ||
        sketch.undo[iundo].var === undefined ||
        sketch.undo[iundo].con === undefined   ) {
        alert("There is nothing to undo");
        return;
    }

    // remove old Sketch info and initialize the arrays
    sketch.pnt        = {};
    sketch.pnt.x      = [];
    sketch.pnt.y      = [];
    sketch.pnt.lbl    = [];

    sketch.seg        = {};
    sketch.seg.type   = [];
    sketch.seg.ibeg   = [];
    sketch.seg.iend   = [];
    sketch.seg.xm     = [];
    sketch.seg.ym     = [];
    sketch.seg.lbl    = [];

    sketch.seg.dip    = [];
    sketch.seg.xc     = [];
    sketch.seg.yc     = [];
    sketch.seg.rad    = [];
    sketch.seg.tbeg   = [];
    sketch.seg.tend   = [];

    sketch.var        = {};
    sketch.var.name   = [];
    sketch.var.value  = [];

    sketch.con        = {};
    sketch.con.type   = [];
    sketch.con.index1 = [];
    sketch.con.index2 = [];
    sketch.con.value  = [];

    // copy the undo information into Sketch information
    sketch.mode = sketch.undo[iundo].mode;

    for (var ipnt = 0; ipnt < sketch.undo[iundo].pnt.x.length; ipnt++) {
        sketch.pnt.x[  ipnt] =  sketch.undo[iundo].pnt.x[  ipnt];
        sketch.pnt.y[  ipnt] =  sketch.undo[iundo].pnt.y[  ipnt];
        sketch.pnt.lbl[ipnt] =  sketch.undo[iundo].pnt.lbl[ipnt];
    }
    for (var iseg = 0; iseg < sketch.undo[iundo].seg.type.length; iseg++) {
        sketch.seg.type[iseg] = sketch.undo[iundo].seg.type[iseg];
        sketch.seg.ibeg[iseg] = sketch.undo[iundo].seg.ibeg[iseg];
        sketch.seg.iend[iseg] = sketch.undo[iundo].seg.iend[iseg];
        sketch.seg.xm[  iseg] = sketch.undo[iundo].seg.xm[  iseg];
        sketch.seg.ym[  iseg] = sketch.undo[iundo].seg.ym[  iseg];
        sketch.seg.lbl[ iseg] = sketch.undo[iundo].seg.lbl[ iseg];
        sketch.seg.dip[ iseg] = sketch.undo[iundo].seg.dip[ iseg];
        sketch.seg.xc[  iseg] = sketch.undo[iundo].seg.xc[  iseg];
        sketch.seg.yc[  iseg] = sketch.undo[iundo].seg.yc[  iseg];
        sketch.seg.rad[ iseg] = sketch.undo[iundo].seg.rad[ iseg];
        sketch.seg.tbeg[iseg] = sketch.undo[iundo].seg.tbeg[iseg];
        sketch.seg.tend[iseg] = sketch.undo[iundo].seg.tend[iseg];
    }
    for (var ivar = 0; ivar <  sketch.undo[iundo].var.name.length; ivar++) {
        sketch.var.name[ ivar] = sketch.undo[iundo].var.name[ ivar];
        sketch.var.value[ivar] = sketch.undo[iundo].var.value[ivar];
    }
    for (var icon = 0; icon <   sketch.undo[iundo].con.type.length; icon++){
        sketch.con.type[  icon] = sketch.undo[iundo].con.type[  icon];
        sketch.con.index1[icon] = sketch.undo[iundo].con.index1[icon];
        sketch.con.index2[icon] = sketch.undo[iundo].con.index2[icon];
        sketch.con.value[ icon] = sketch.undo[iundo].con.value[ icon];
    }

    // remove the undo information
    sketch.undo[iundo] = {};

    // decrement number of undos
    sketch.nundo--;

    // update the display
    sketcherDraw();
};


//
// save an undo snapshot
//
var sketcherSaveUndo = function () {
    // alert("sketcherSaveUndo()");

    var iundo = sketch.nundo % sketch.maxundo;

    // remove old undo info and initialize the arrays
    sketch.undo[iundo]            = {};

    sketch.undo[iundo].pnt        = {};
    sketch.undo[iundo].pnt.x      = [];
    sketch.undo[iundo].pnt.y      = [];
    sketch.undo[iundo].pnt.lbl    = [];

    sketch.undo[iundo].seg        = {};
    sketch.undo[iundo].seg.type   = [];
    sketch.undo[iundo].seg.ibeg   = [];
    sketch.undo[iundo].seg.iend   = [];
    sketch.undo[iundo].seg.xm     = [];
    sketch.undo[iundo].seg.ym     = [];
    sketch.undo[iundo].seg.lbl    = [];

    sketch.undo[iundo].seg.dip    = [];
    sketch.undo[iundo].seg.xc     = [];
    sketch.undo[iundo].seg.yc     = [];
    sketch.undo[iundo].seg.rad    = [];
    sketch.undo[iundo].seg.tbeg   = [];
    sketch.undo[iundo].seg.tend   = [];

    sketch.undo[iundo].var        = {};
    sketch.undo[iundo].var.name   = [];
    sketch.undo[iundo].var.value  = [];

    sketch.undo[iundo].con        = {};
    sketch.undo[iundo].con.type   = [];
    sketch.undo[iundo].con.index1 = [];
    sketch.undo[iundo].con.index2 = [];
    sketch.undo[iundo].con.value  = [];

    // copy the Sketch information into the undo storage
    sketch.undo[iundo].mode = sketch.mode;

    for (var ipnt = 0; ipnt < sketch.pnt.x.length; ipnt++) {
        sketch.undo[iundo].pnt.x[  ipnt] = sketch.pnt.x[  ipnt];
        sketch.undo[iundo].pnt.y[  ipnt] = sketch.pnt.y[  ipnt];
        sketch.undo[iundo].pnt.lbl[ipnt] = sketch.pnt.lbl[ipnt];
    }
    for (var iseg = 0; iseg < sketch.seg.type.length; iseg++) {
        sketch.undo[iundo].seg.type[iseg] = sketch.seg.type[iseg];
        sketch.undo[iundo].seg.ibeg[iseg] = sketch.seg.ibeg[iseg];
        sketch.undo[iundo].seg.iend[iseg] = sketch.seg.iend[iseg];
        sketch.undo[iundo].seg.xm[  iseg] = sketch.seg.xm[  iseg];
        sketch.undo[iundo].seg.ym[  iseg] = sketch.seg.ym[  iseg];
        sketch.undo[iundo].seg.lbl[ iseg] = sketch.seg.lbl[ iseg];
        sketch.undo[iundo].seg.dip[ iseg] = sketch.seg.dip[ iseg];
        sketch.undo[iundo].seg.xc[  iseg] = sketch.seg.xc[  iseg];
        sketch.undo[iundo].seg.yc[  iseg] = sketch.seg.yc[  iseg];
        sketch.undo[iundo].seg.rad[ iseg] = sketch.seg.rad[ iseg];
        sketch.undo[iundo].seg.tbeg[iseg] = sketch.seg.tbeg[iseg];
        sketch.undo[iundo].seg.tend[iseg] = sketch.seg.tend[iseg];
    }
    for (var ivar = 0; ivar < sketch.var.name.length; ivar++) {
        sketch.undo[iundo].var.name[ ivar] = sketch.var.name[ ivar];
        sketch.undo[iundo].var.value[ivar] = sketch.var.value[ivar];
    }
    for (var icon = 0; icon < sketch.con.type.length; icon++){
        sketch.undo[iundo].con.type[  icon] = sketch.con.type[  icon];
        sketch.undo[iundo].con.index1[icon] = sketch.con.index1[icon];
        sketch.undo[iundo].con.index2[icon] = sketch.con.index2[icon];
        sketch.undo[iundo].con.value[ icon] = sketch.con.value[ icon];
    }

    // increment number of undos
    sketch.nundo++;
};


//
// load a Sketch from a semicolon-separated strings sent from serveCSM
//
var sketcherLoad = function (begs, vars, cons, segs) {
    // alert("sketcherLoad()");
    // alert("   begs="+begs);
    // alert("   vars="+vars);
    // alert("   cons="+cons);
    // alert("   segs="+segs);

    // break up the inputs
    var begsList = begs.split(";");
    var varsList = vars.split(";");
    var consList = cons.split(";");
    var segsList = segs.split(";");

    // get the number of points, segments, variables, and constraints
    var npnt = (varsList.length - 1) / 3;
    var ncon = (consList.length - 1) / 4;
    var nseg = (segsList.length - 1) / 3;

    // store the base expressions from the SKBEG statement
    sketch.skbegX   = begsList[0];
    sketch.skbegY   = begsList[2];
    sketch.skbegZ   = begsList[4];
    sketch.relative = begsList[6];
    sketch.suggest  = undefined;

    // find the extrema of the coordinates
    var xmin = Number(begsList[1]);
    var xmax = Number(begsList[1]);
    var ymin = Number(begsList[3]);
    var ymax = Number(begsList[3]);

    for (var ipnt = 1; ipnt < npnt; ipnt++) {
        xmin = Math.min(xmin, Number(varsList[3*ipnt  ]));
        xmax = Math.max(xmax, Number(varsList[3*ipnt  ]));
        ymin = Math.min(ymin, Number(varsList[3*ipnt+1]));
        ymax = Math.max(ymax, Number(varsList[3*ipnt+1]));
    }

    // get the size of the canvas
    var canvas = document.getElementById("sketcher");
    var width  = canvas.clientWidth;
    var height = canvas.clientHeight;

    if (npnt > 1) {
        var test1 = (xmax - xmin) / width;
        var test2 = (ymax - ymin) / height;

        // set up sizing so that Sketch fills middle 50 percent of the window
        sketch.scale = 2 * Math.max(test1, test2);
        sketch.xorig = (width  - (xmin + xmax) / sketch.scale) / 2;
        sketch.yorig = (height + (ymin + ymax) / sketch.scale) / 2;

    // not enough data to set scales, so guess something base upon the SKBEG statement
    } else {
        sketch.scale = 0.01;
        sketch.xorig = width  / 2 - begsList[1] / sketch.scale;
        sketch.yorig = height / 2 + begsList[3] / sketch.scale;
    }

    // create the points
    sketch.pnt.x   = [];
    sketch.pnt.y   = [];
    sketch.pnt.lbl = [];

    for (ipnt = 0; ipnt < npnt; ipnt++) {
//        sketch.pnt.x[  ipnt] = Math.floor(sketch.xorig + Number(varsList[3*ipnt  ]) / sketch.scale);
//        sketch.pnt.y[  ipnt] = Math.floor(sketch.yorig - Number(varsList[3*ipnt+1]) / sketch.scale);
        sketch.pnt.x[  ipnt] = sketch.xorig + Number(varsList[3*ipnt  ]) / sketch.scale;
        sketch.pnt.y[  ipnt] = sketch.yorig - Number(varsList[3*ipnt+1]) / sketch.scale;
        sketch.pnt.lbl[ipnt] = "";
    }

    // overwrite the first point by the begsList
    if (npnt > 0) {
//        sketch.pnt.x[0] = Math.floor(sketch.xorig + Number(begsList[1]) / sketch.scale);
//        sketch.pnt.y[0] = Math.floor(sketch.yorig - Number(begsList[3]) / sketch.scale);
        sketch.pnt.x[0] = sketch.xorig + Number(begsList[1]) / sketch.scale;
        sketch.pnt.y[0] = sketch.yorig - Number(begsList[3]) / sketch.scale;
    }

    // create the segments
    sketch.seg.type = [];
    sketch.seg.ibeg = [];
    sketch.seg.iend = [];
    sketch.seg.xm   = [];
    sketch.seg.ym   = [];
    sketch.seg.lbl  = [];

    sketch.seg.dip  = [];
    sketch.seg.xc   = [];
    sketch.seg.yc   = [];
    sketch.seg.rad  = [];
    sketch.seg.tbeg = [];
    sketch.seg.tend = [];

    for (var iseg = 0; iseg < nseg; iseg++) {
        sketch.seg.type[iseg] =        segsList[3*iseg  ];
        sketch.seg.ibeg[iseg] = Number(segsList[3*iseg+1]) - 1;
        sketch.seg.iend[iseg] = Number(segsList[3*iseg+2]) - 1;
        sketch.seg.lbl[ iseg] = "";

        sketch.seg.dip[ iseg] = 0;
        sketch.seg.xc[  iseg] = 0;
        sketch.seg.yc[  iseg] = 0;
        sketch.seg.rad[ iseg] = 0;
        sketch.seg.tbeg[iseg] = 0;
        sketch.seg.tend[iseg] = 0;

        // update the segment if an arc
        if (sketch.seg.type[iseg] == "C") {
            if (iseg < npnt-1) {
                sketch.seg.dip[iseg] = Number(varsList[3*iseg+5]) / sketch.scale;
            } else {
                sketch.seg.dip[iseg] = Number(varsList[       2]) / sketch.scale;
            }
        }

        sketcherUpdateSegData(iseg);
    }

    // create the variables
    sketch.var.name  = [];
    sketch.var.value = [];

    var nvar = 0;
    for (var ivar = 0; ivar < npnt; ivar++) {
        sketch.var.name[ nvar] = sprintf("::x[%]", ivar+1);
        sketch.var.value[nvar] = varsList[3*ivar  ];
        nvar++;

        sketch.var.name[ nvar] = sprintf("::y[%]", ivar+1);
        sketch.var.value[nvar] = varsList[3*ivar+1];
        nvar++;

        if (varsList.length > 0 && Number(varsList[3*ivar+2]) != 0) {
            sketch.var.name[ nvar] = sprintf("::d[%]", ivar+1);
            sketch.var.value[nvar] = varsList[3*ivar+2];
            nvar++;
        }
    }

    // create the constraints
    sketch.con.type   = [];
    sketch.con.index1 = [];
    sketch.con.index2 = [];
    sketch.con.value  = [];

    for (var icon = 0; icon < ncon; icon++) {
        sketch.con.type[  icon] =        consList[4*icon  ];
        sketch.con.index1[icon] = Number(consList[4*icon+1]) - 1;
        sketch.con.index2[icon] = Number(consList[4*icon+2]) - 1;
        sketch.con.value[ icon] =        consList[4*icon+3];

        if (sketch.con.type[icon] == "W" || sketch.con.type[icon] == "D") {
            // not associated with Point or Segment
        } else if (sketch.con.index2[icon] < 0) {
            sketch.con.index2[icon] += 1;
            sketch.pnt.lbl[sketch.con.index1[icon]] += sketch.con.type[icon];
        } else {
            sketch.seg.lbl[sketch.con.index1[icon]] += sketch.con.type[icon];
        }
    }

    // if there were no constraints, create X and Y constraints at first Point
    if (ncon == 0) {
//        sketch.pnt.x[  0] = Math.floor(sketch.xorig + Number(begsList[1]) / sketch.scale);
//        sketch.pnt.y[  0] = Math.floor(sketch.yorig - Number(begsList[3]) / sketch.scale);
        sketch.pnt.x[  0] = sketch.xorig + Number(begsList[1]) / sketch.scale;
        sketch.pnt.y[  0] = sketch.yorig - Number(begsList[3]) / sketch.scale;
        sketch.pnt.lbl[0] = "XY";

        sketch.var.name[ 0] = "::x[1]";
        sketch.var.value[0] = sketch.pnt.x[0];

        sketch.var.name[ 1] = "::y[1]";
        sketch.var.value[1] = sketch.pnt.y[0];

        if (sketch.relative == 0) {
            sketch.con.type[  0] = "X";
            sketch.con.index1[0] = 0;
            sketch.con.index2[0] = -1;
            sketch.con.value[ 0] = sketch.skbegX;

            sketch.con.type[  1] = "Y";
            sketch.con.index1[1] = 0;
            sketch.con.index2[1] = -1;
            sketch.con.value[ 1] = sketch.skbegY;
        } else {
            sketch.con.type[  0] = "X";
            sketch.con.index1[0] = 0;
            sketch.con.index2[0] = -1;
            sketch.con.value[ 0] = 0;

            sketch.con.type[  1] = "Y";
            sketch.con.index1[1] = 0;
            sketch.con.index2[1] = -1;
            sketch.con.value[ 1] = 0;
        }
    }

    // set the mode depending on status of Sketch
    if (sketch.seg.type.length == 0) {
        sketch.mode = 1;
    } else if (sketch.seg.iend[sketch.seg.iend.length-1] != 0) {
        sketch.mode = 1;
    } else {
        sketch.mode = 3;
    }

    // draw the current Sketch
    sketcherDraw();
};


//
// update the variables associated with a segment
//
var sketcherUpdateSegData = function (iseg) {
    // alert("sketcherUpdateSegData(iseg="+iseg+")");

    sketch.seg.xm[iseg] = (sketch.pnt.x[sketch.seg.ibeg[iseg]] + sketch.pnt.x[sketch.seg.iend[iseg]]) / 2;
    sketch.seg.ym[iseg] = (sketch.pnt.y[sketch.seg.ibeg[iseg]] + sketch.pnt.y[sketch.seg.iend[iseg]]) / 2;

    if (Math.abs(sketch.seg.dip[iseg]) > 1.0e-6) {
        var xa = sketch.pnt.x[sketch.seg.ibeg[iseg]];
        var ya = sketch.pnt.y[sketch.seg.ibeg[iseg]];
        var xb = sketch.pnt.x[sketch.seg.iend[iseg]];
        var yb = sketch.pnt.y[sketch.seg.iend[iseg]];
        var d  =            sketch.seg.dip[ iseg];

        var L  = Math.sqrt((xb-xa) * (xb-xa) + (yb-ya) * (yb-ya));
        var R  = (L*L + 4*d*d) / (8*d);

        sketch.seg.xc[iseg] = (xa + xb) / 2 + (R - d) * (yb - ya) / L;
        sketch.seg.yc[iseg] = (ya + yb) / 2 - (R - d) * (xb - xa) / L;

        sketch.seg.rad[ iseg] = R;
        sketch.seg.tbeg[iseg] = Math.atan2(ya-sketch.seg.yc[iseg], xa-sketch.seg.xc[iseg]);
        sketch.seg.tend[iseg] = Math.atan2(yb-sketch.seg.yc[iseg], xb-sketch.seg.xc[iseg]);

        if (R > 0) {
            sketch.seg.tbeg[iseg] -= 2 * Math.PI;
            while (sketch.seg.tbeg[iseg] < sketch.seg.tend[iseg]) {
                sketch.seg.tbeg[iseg] += 2 * Math.PI;
            }
        } else {
            sketch.seg.tbeg[iseg] += 2 * Math.PI;
            while (sketch.seg.tbeg[iseg] > sketch.seg.tend[iseg]) {
                sketch.seg.tbeg[iseg] -= 2 * Math.PI;
            }
        }

        var tavg = (sketch.seg.tbeg[iseg] + sketch.seg.tend[iseg]) / 2;
        sketch.seg.xm[iseg] = sketch.seg.xc[iseg] + Math.cos(tavg) * Math.abs(R);
        sketch.seg.ym[iseg] = sketch.seg.yc[iseg] + Math.sin(tavg) * Math.abs(R);
    }
};


//
// Draw the current Sketch
//
var sketcherDraw = function () {
    // alert("sketcherDraw()");

    // set up context
    var canvas  = document.getElementById("sketcher");
    var context = canvas.getContext("2d");

    // clear the screen
    context.clearRect(0, 0, canvas.width, canvas.height);

    context.lineCap  = "round";
    context.lineJoin = "round";

    // determine number of Points and Segments in Sketch
    var npnt = sketch.pnt.x.length;
    var nseg = sketch.seg.type.length;
    var ncon = sketch.con.type.length;

    // return if no Points
    if (npnt <= 0) {
        var button = document.getElementById("solveButton");
        button.style.backgroundColor =  "#FFFF3F";
        button["innerHTML"] = "Initializing...";

        return;
    }

    // axes
    context.strokeStyle = "yellow";
    context.lineWidth   = 1;

    context.beginPath();
    context.moveTo(sketch.pnt.x[0], -100);
    context.lineTo(sketch.pnt.x[0], 1000);
    context.stroke();

    context.beginPath();
    context.moveTo(-100, sketch.pnt.y[0]);
    context.lineTo(1000, sketch.pnt.y[0]);
    context.stroke();

    // draw the part of the Sketch that is complete
    if (nseg > 0) {

        // we are editting the Cirarc at the end, so it has to
        //    be drawn separately (since it will be drawn in red)
        if (sketch.mode == 2) {

            // the part of the Sketch that is done
            context.strokeStyle = "black";
            context.lineWidth   = 3;
            context.beginPath();

            var ibeg = sketch.seg.ibeg[0];
            var iend = sketch.seg.iend[0];
            context.moveTo(sketch.pnt.x[ibeg], sketch.pnt.y[ibeg]);

            for (var iseg = 0; iseg < nseg-1; iseg++) {
                iend = sketch.seg.iend[iseg];

                if (Math.abs(sketch.seg.dip[iseg]) < 1e-6) {
                    context.lineTo(sketch.pnt.x[iend], sketch.pnt.y[iend]);
                } else if (sketch.seg.rad[iseg] > 0) {
                    context.arc(sketch.seg.xc[iseg], sketch.seg.yc[iseg],
                                +sketch.seg.rad[iseg], sketch.seg.tbeg[iseg], sketch.seg.tend[iseg], 1);
                } else {
                    context.arc(sketch.seg.xc[iseg], sketch.seg.yc[iseg],
                                -sketch.seg.rad[iseg], sketch.seg.tbeg[iseg], sketch.seg.tend[iseg], 0);
                }
            }

            context.stroke();

            // the Cirarc at the end which is to be rendered in red
            context.strokeStyle = "red";
            context.beginPath();

            if (sketch.seg.dip[nseg-1] > 0) {
                context.arc(sketch.seg.xc[nseg-1], sketch.seg.yc[nseg-1],
                            +sketch.seg.rad[nseg-1], sketch.seg.tbeg[nseg-1], sketch.seg.tend[nseg-1], 1);
            } else {
                context.arc(sketch.seg.xc[nseg-1], sketch.seg.yc[nseg-1],
                            -sketch.seg.rad[nseg-1], sketch.seg.tbeg[nseg-1], sketch.seg.tend[nseg-1], 0);
            }

            context.stroke();

        // we are NOT editting the Cirarc at the end
        } else {
            context.strokeStyle = "black";
            context.lineWidth   = 3;
            context.beginPath();

            var ibeg = sketch.seg.ibeg[0];
            var iend = sketch.seg.iend[0];
            context.moveTo(sketch.pnt.x[ibeg], sketch.pnt.y[ibeg]);

            for (var iseg = 0; iseg < nseg; iseg++) {
                iend = sketch.seg.iend[iseg];

                if (sketch.seg.type[iseg] == "B") {
                    var Xbezier = [];
                    var Ybezier = [];
                    var nbezier = 2;

                    Xbezier[0] = sketch.pnt.x[ibeg];
                    Ybezier[0] = sketch.pnt.y[ibeg];
                    Xbezier[1] = sketch.pnt.x[iend];
                    Ybezier[1] = sketch.pnt.y[iend];

                    while (iseg < nseg-1) {
                        if (sketch.seg.type[iseg+1] == "B") {
                            iseg++;
                            iend = sketch.seg.iend[iseg];
                            Xbezier[nbezier] = sketch.pnt.x[iend];
                            Ybezier[nbezier] = sketch.pnt.y[iend];
                            nbezier++;
                        } else {
                            break;
                        }
                    }

                    sketcherDrawBezier(context, Xbezier, Ybezier);
//$$$             } else if (sketch.seg.type[iseg] == "S") {
//$$$                 var Xspline = [];
//$$$                 var Yspline = [];
//$$$                 var nspline = 2;
//
//$$$                 Xspline[0] = sketch.pnt.x[ibeg];
//$$$                 Yspline[0] = sketch.pnt.y[ibeg];
//$$$                 Xspline[1] = sketch.pnt.x[iend];
//$$$                 Yspline[1] = sketch.pnt.y[iend];
//
//$$$                 while (iseg < nseg-1) {
//$$$                     if (sketch.seg.type[iseg+1] == "S") {
//$$$                         iseg++;
//$$$                         iend = sketch.seg.iend[iseg];
//$$$                         Xspline[nspline] = sketch.pnt.x[iend];
//$$$                         Yspline[nspline] = sketch.pnt.y[iend];
//$$$                         nspline++;
//$$$                     } else {
//$$$                         break;
//$$$                     }
//$$$                 }
//
//$$$                 sketcherDrawSpline(context, Xspline, Yspline);
                } else if (Math.abs(sketch.seg.dip[iseg]) < 1e-6) {
                    context.lineTo(sketch.pnt.x[iend], sketch.pnt.y[iend]);
                } else if (sketch.seg.dip[iseg] > 0) {
                    context.arc(sketch.seg.xc[iseg], sketch.seg.yc[iseg],
                                +sketch.seg.rad[iseg], sketch.seg.tbeg[iseg], sketch.seg.tend[iseg], 1);
                } else {
                    context.arc(sketch.seg.xc[iseg], sketch.seg.yc[iseg],
                                -sketch.seg.rad[iseg], sketch.seg.tbeg[iseg], sketch.seg.tend[iseg], 0);
                }

                ibeg = iend;
            }

            context.stroke();

            // if Sketch is closed, fill it
            if (iend == sketch.seg.ibeg[0]) {
                if (sketch.var.name.length == sketch.con.type.length) {
                    context.fillStyle = "#c0f0c0";  // light green
                    context.fill();
                } else {
                    context.fillStyle = "#e0e0e0";  // grey
                    context.fill();
                }
            }

            // if draw straight lines for Segments associated with splines
            //    and beziers
            context.strokeStyle = "cyan";
            context.lineWidth   = 1;
            for (iseg = 0; iseg < nseg; iseg++) {
                if (sketch.seg.type[iseg] == "B" ||
                    sketch.seg.type[iseg] == "S"   ) {

                    ibeg = sketch.seg.ibeg[iseg];
                    iend = sketch.seg.iend[iseg];
                    context.beginPath();
                    context.moveTo(sketch.pnt.x[ibeg], sketch.pnt.y[ibeg]);
                    context.lineTo(sketch.pnt.x[iend], sketch.pnt.y[iend]);
                    context.stroke();
                }
            }

            // we are creating W constraint
            if (sketch.mode == 4) {
                var xmid = (  sketch.pnt.x[sketch.basePoint] + wv.cursorX) / 2;
                var ymid = (2*sketch.pnt.y[sketch.basePoint] + wv.cursorY) / 3;

                context.strokeStyle = "red";
                context.lineWidth   = 1;
                context.beginPath();
                context.moveTo(sketch.pnt.x[sketch.basePoint], sketch.pnt.y[sketch.basePoint]);
                context.lineTo(sketch.pnt.x[sketch.basePoint], ymid                      );
                context.lineTo(wv.cursorX,                 ymid                      );
                context.lineTo(wv.cursorX,                 wv.cursorY                );
                context.stroke();
            }

            // we are creating D constraint
            if (sketch.mode == 5) {
                var xmid = (2*sketch.pnt.x[sketch.basePoint] + wv.cursorX) / 3;
                var ymid = (  sketch.pnt.y[sketch.basePoint] + wv.cursorY) / 2;

                context.strokeStyle = "red";
                context.lineWidth   = 1;
                context.beginPath();
                context.moveTo(sketch.pnt.x[sketch.basePoint], sketch.pnt.y[sketch.basePoint]);
                context.lineTo(xmid,                       sketch.pnt.y[sketch.basePoint]);
                context.lineTo(xmid,                       wv.cursorY                );
                context.lineTo(wv.cursorX,                 wv.cursorY                );
                context.stroke();
            }
        }
    }

    // draw the Points
    context.fillStyle = "black";
    context.fillRect(sketch.pnt.x[0]-3, sketch.pnt.y[0]-3, 7, 7);

    for (var ipnt = 1; ipnt < npnt; ipnt++) {
        context.fillRect(sketch.pnt.x[ipnt]-2, sketch.pnt.y[ipnt]-2, 5, 5);
    }

    // display constraint labels associated with the Points
    for (ipnt = 0; ipnt < npnt; ipnt++) {
        if (sketch.pnt.lbl[ipnt] != "") {
            context.font         = "10px Verdana";
            context.textAlign    = "center";
            context.textBaseline = "middle";

            var width  = context.measureText(sketch.pnt.lbl[ipnt]).width;
            var height = 8;

            context.fillStyle    = "yellow";
            context.fillRect(sketch.pnt.x[ipnt]-width/2-2, sketch.pnt.y[ipnt]-height/2-2, width+4, height+4);

            context.fillStyle    = "black";
            context.fillText(sketch.pnt.lbl[ipnt], sketch.pnt.x[ipnt], sketch.pnt.y[ipnt]);
        }
    }

    // display constraint labels associated with the Segments
    for (var iseg = 0; iseg < nseg; iseg++) {
        if (sketch.seg.lbl[iseg] != "") {
            context.font         = "12px Verdana";
            context.textAlign    = "center";
            context.textBaseline = "middle";

            var width  = context.measureText(sketch.seg.lbl[iseg]).width;
            var height = 8;

            context.fillStyle    = "yellow";
            context.fillRect(sketch.seg.xm[iseg]-width/2-2, sketch.seg.ym[iseg]-height/2-2, width+4, height+4);

            context.fillStyle    = "black";
            context.fillText(sketch.seg.lbl[iseg], sketch.seg.xm[iseg], sketch.seg.ym[iseg]);
        }
    }

    // display width and depth constraints
    for (var icon = 0; icon < ncon; icon++) {
        if (sketch.con.type[icon] == "W" || sketch.con.type[icon] == "D") {
            var ibeg = sketch.con.index1[icon];
            var iend = sketch.con.index2[icon];

            context.strokeStyle = "green";
            context.lineWidth   = 1;
            context.beginPath();
            context.moveTo(sketch.pnt.x[ibeg], sketch.pnt.y[ibeg]);
            if (sketch.con.type[icon] == "W") {
                var xmid = (  sketch.pnt.x[ibeg] + sketch.pnt.x[iend]) / 2;
                var ymid = (2*sketch.pnt.y[ibeg] + sketch.pnt.y[iend]) / 3;

                context.lineTo(sketch.pnt.x[ibeg], ymid);
                context.lineTo(sketch.pnt.x[iend], ymid);
            } else {
                var xmid = (2*sketch.pnt.x[ibeg] + sketch.pnt.x[iend]) / 3;
                var ymid = (  sketch.pnt.y[ibeg] + sketch.pnt.y[iend]) / 2;

                context.lineTo(xmid, sketch.pnt.y[ibeg]);
                context.lineTo(xmid, sketch.pnt.y[iend]);
            }
            context.lineTo(sketch.pnt.x[iend], sketch.pnt.y[iend]);
            context.stroke();

            context.font         = "12px Verdana";
            context.textAlign    = "center";
            context.textBaseline = "middle";

            if (sketch.con.type[icon] == "W") {
                var label = "W";
            } else {
                var label = "D";
            }
            var width  = context.measureText(label).width;
            var height = 8;

            context.fillStyle = "yellow";
            context.fillRect(xmid-width/2-2, ymid-height/2-2, width+4, height+4);

            context.fillStyle = "black";
            context.fillText(label, xmid, ymid);
        }
    }

    if (sketch.suggest) {
        // suggested Sketch deletions
        if        (sketch.suggest.substring(0,4) == "*del") {
            var suggestions = sketch.suggest.split(";");

            for (ipnt = 0; ipnt < sketch.pnt.x.length; ipnt++) {
                var lbl = "";

                for (var ient = 1; ient < suggestions.length; ient+=3) {
                    if (suggestions[ient+1] == ipnt+1 &&
                        suggestions[ient+2] == -1     ) {
                        lbl += suggestions[ient];
                    }
                }

                if (lbl.length > 0) {
                    context.font         = "10px Verdana";
                    context.textAlign    = "center";
                    context.textBaseline = "middle";

                    var width  = context.measureText(lbl).width;
                    var height = 8;

                    context.fillStyle    = "red";

                    context.clearRect(sketch.pnt.x[ipnt]-width/2-2, sketch.pnt.y[ipnt]-2*height-2, width+4, height+4);
                    context.fillText(lbl, sketch.pnt.x[ipnt], sketch.pnt.y[ipnt]-3*height/2);
                }
            }

            for (ipnt = 0; ipnt < sketch.pnt.x.length; ipnt++) {
                var lbl = "";

                for (var ient = 1; ient < suggestions.length; ient+=3) {
                    if (suggestions[ient+1] == ipnt+1 &&
                        suggestions[ient+2] >  -1       ) {
                        lbl += suggestions[ient];
                    }
                }

                if (lbl.length > 0) {
                    context.font         = "10px Verdana";
                    context.textAlign    = "center";
                    context.textBaseline = "middle";

                    var width  = context.measureText(lbl).width;
                    var height = 8;

                    context.fillStyle    = "red";

                    context.clearRect(sketch.seg.xm[ipnt]-width/2-2, sketch.seg.ym[ipnt]-2*height-2, width+4, height+4);
                    context.fillText(lbl, sketch.seg.xm[ipnt], sketch.seg.ym[ipnt]-3*height/2);
                }
            }

            return;

        // suggested Sketch additions
        } else if (sketch.suggest.substring(0,4) == "*add") {
            var suggestions = sketch.suggest.split(";");

            for (ipnt = 0; ipnt < sketch.pnt.x.length; ipnt++) {
                var lbl = "";

                for (var ient = 1; ient < suggestions.length; ient+=3) {
                    if (suggestions[ient+1] == ipnt+1 &&
                        suggestions[ient+2] == -1     ) {
                        lbl += suggestions[ient];
                    }
                }

                if (lbl.length > 0) {
                    context.font         = "10px Verdana";
                    context.textAlign    = "center";
                    context.textBaseline = "middle";

                    var width  = context.measureText(lbl).width;
                    var height = 8;

                    context.fillStyle    = "green";

                    context.clearRect(sketch.pnt.x[ipnt]-width/2-2, sketch.pnt.y[ipnt]-2*height-2, width+4, height+4);
                    context.fillText(lbl, sketch.pnt.x[ipnt], sketch.pnt.y[ipnt]-3*height/2);
                }
            }

            for (ipnt = 0; ipnt < sketch.pnt.x.length; ipnt++) {
                var lbl = "";

                for (var ient = 1; ient < suggestions.length; ient+=3) {
                    if (suggestions[ient+1] == ipnt+1 &&
                        suggestions[ient+2] >  -1       ) {
                        lbl += suggestions[ient];
                    }
                }

                if (lbl.length > 0) {
                    context.font         = "10px Verdana";
                    context.textAlign    = "center";
                    context.textBaseline = "middle";

                    var width  = context.measureText(lbl).width;
                    var height = 8;

                    context.fillStyle    = "green";

                    context.clearRect(sketch.seg.xm[ipnt]-width/2-2, sketch.seg.ym[ipnt]-2*height-2, width+4, height+4);
                    context.fillText(lbl, sketch.seg.xm[ipnt], sketch.seg.ym[ipnt]-3*height/2);
                }
            }

            return;
        }
    }

    // if we are in mode 1, draw a line from the last Point to the cursor
    if (sketch.mode == 1) {

        // nominal line to cursor is blue
        var color = "blue";

        // if horizontal, make it orange
        if (Math.abs(sketch.pnt.y[npnt-1]-wv.cursorY) < sketch.halo) {
            color = "orange";

        // otherwise see if it alignes with another point
        } else {
            for (var jpnt = 0; jpnt < npnt-1; jpnt++) {
                if (Math.abs(sketch.pnt.y[jpnt]-wv.cursorY) < sketch.halo) {
                    context.strokeStyle = "orange";
                    context.lineWidth   = 1;

                    context.beginPath();
                    context.moveTo(sketch.pnt.x[jpnt], sketch.pnt.y[jpnt]);
                    context.lineTo(wv.cursorX,  wv.cursorY );
                    context.stroke();
                    break;
                }
            }
        }

        // if vertical, make it orange
        if (Math.abs(sketch.pnt.x[npnt-1]-wv.cursorX) < sketch.halo) {
            color = "orange";

        // otherwise see if it alignes with another point
        } else {
            for (var jpnt = 0; jpnt < npnt-1; jpnt++) {
                if (Math.abs(sketch.pnt.x[jpnt]-wv.cursorX) < sketch.halo) {
                    context.strokeStyle = "orange";
                    context.lineWidth   = 1;

                    context.beginPath();
                    context.moveTo(sketch.pnt.x[jpnt], sketch.pnt.y[jpnt]);
                    context.lineTo(wv.cursorX,  wv.cursorY );
                    context.stroke();
                    break;
                }
            }
        }

        // if the cursor is close to the first Point, indicate that
        if (Math.abs(sketch.pnt.x[0]-wv.cursorX) < sketch.halo &&
            Math.abs(sketch.pnt.y[0]-wv.cursorY) < sketch.halo   ) {
            context.strokeStyle = "orange";
            context.lineWidth   = 5;
            context.beginPath();
            context.arc(sketch.pnt.x[0], sketch.pnt.y[0], 17, 0, 2*Math.PI);
            context.stroke();
        }

        // draw the line to the cursor
        context.strokeStyle = color;
        context.lineWidth   = 2;

        context.beginPath();
        context.moveTo(sketch.pnt.x[npnt-1], sketch.pnt.y[npnt-1]);
        context.lineTo(wv.cursorX,    wv.cursorY   );
        context.stroke();
    }

    // have the Build button show the current Sketcher status
    var button = document.getElementById("solveButton");
    button.style.backgroundColor =  "#FFFF3F";

    if        (sketch.mode == 0) {
        button["innerHTML"] = "Initializing...";
    } else if (sketch.mode == 1) {
        button["innerHTML"] = "Drawing...";
    } else if (sketch.mode == 2) {
        button["innerHTML"] = "Setting R...";
//  } else if (sketch.mode == 3) {
//      handled below
    } else if (sketch.mode == 4) {
        button["innerHTML"] = "Setting W...";
    } else if (sketch.mode == 5) {
        button["innerHTML"] = "Setting D...";
    } else if (sketch.mode == 6) {
        button["innerHTML"] = "Up to date";
        button.style.backgroundColor = null;
        document.getElementById("doneMenuBtn").style.backgroundColor = "#3FFF3F";
    } else if (sketch.var.name.length != sketch.con.type.length) {
        button["innerHTML"] = "Constraining...";
        button.style.backgroundColor = "#FFFF3F"
        document.getElementById("doneMenuBtn").style.backgroundColor = null;
    } else {
        button["innerHTML"] = "Press to Solve";
        button.style.backgroundColor = "#3FFF3F";
    }

    // post informtion about current mode in blframe
    var skstat = document.getElementById("sketcherStatus");

    var mesg = "ndof=" + sketch.var.name.length + "   ncon=" + sketch.con.type.length + "\n";
    if        (sketch.mode == 1) {
        mesg += "Valid commands are:\n";
        mesg += "  'l'   add linseg\n";
        mesg += "  'c'   add cirarc\n";
        mesg += "  's'   add spline\n";
        mesg += "  'b'   add bezier\n";
        mesg += "  'z'   add zero-length segment\n";
        mesg += "  'o'   complete (open) sketch\n";
    } else if (sketch.mode == 2) {
        mesg += "Hit any character to set curvature\n";
    } else if (sketch.mode == 3) {
        mesg += "Valid constraints at points\n";
        mesg += "  'x' (fix x)     'y' (fix y)\n";
        mesg += "  'p' (perp)      't' (tangent)\n";
        mesg += "  'a' (angle)\n";
        mesg += "  'w' (width)     'd' (depth)\n";
        mesg += "Valid constraints on segments\n";
        mesg += "  'h' (horiz)     'v' (vertical)\n";
        mesg += "  'i' (incline)   'l' (length)\n";
        mesg += "Valid constraints on cirarcs\n";
        mesg += "  'r' (radius)    's' (sweep angle)\n";
        mesg += "  'x' (x center)  'y' (y center)\n";
        mesg += "Valid commands anywhere\n";
        mesg += "  '<' (delete)    '?' (info)\n";
    } else if (sketch.mode == 4) {
        mesg += "Hit any character to set width\n";
    } else if (sketch.mode == 5) {
        mesg += "Hit any character to set depth\n";
    }

    var pre  = document.createElement("pre");
    var text = document.createTextNode(mesg);
    pre.appendChild(text);

    skstat.replaceChild(pre, skstat.lastChild);
};


// draw a Bezier curve (with evaluations via Hoerner-like scheme)
//
var sketcherDrawBezier = function (context, Xbezier, Ybezier) {
    // alert("sketcherDrawBezier()");

    var degree = Xbezier.length - 1;

    for (var j = 1; j <= 4*degree; j++) {
        var t  = j / (4.0 * degree);
        var t1 = 1.0 - t;

        var fact   = 1.0;
        var bicoef = 1;

        var xx = Xbezier[0] * t1;
        var yy = Ybezier[0] * t1;

        for (var i = 1; i < degree; i++) {
            fact    *= t;
            bicoef  *= (degree - i + 1) / i;

            xx = (xx + fact * bicoef * Xbezier[i]) * t1;
            yy = (yy + fact * bicoef * Ybezier[i]) * t1;
        }

        xx += fact * t * Xbezier[degree];
        yy += fact * t * Ybezier[degree];

        context.lineTo(xx, yy);
    }
};


// draw a Spline curve (with evaluations via Hoerner-like scheme)
//
var sketcherDrawSpline = function (context, Xspline, Yspline) {
    // alert("sketcherDrawSpline()");

    var degree = Xspline.length - 1;

    // convert Spline points into cubic Bezier control points
    var Xbezier = [];
    var Ybezier = [];
    var gam     = [];

    var bet    = 4;
    Xbezier[1] = (6 * Xspline[1] - Xspline[0]) / bet;
    Ybezier[1] = (6 * Yspline[1] - Yspline[0]) / bet;

    for (var i = 1; i < degree-2; i++) {
        gam[i]       = 1 / bet;
        bet          =  4 - gam[i];
        Xbezier[i+1] = (6 * Xspline[i+1] - Xbezier[i]) / bet;
        Ybezier[i+1] = (6 * Yspline[i+1] - Ybezier[i]) / bet;
    }

    gam[degree-2]     = 1 / bet;
    bet               =  4 - gam[degree-2];
    Xbezier[degree-1] = (6 * Xspline[degree-1] - Xspline[degree] - Xbezier[degree-2]) / bet;
    Ybezier[degree-1] = (6 * Yspline[degree-1] - Yspline[degree] - Ybezier[degree-2]) / bet;

    for (i = degree-2; i > 0; i--) {
        Xbezier[i] -= gam[i] * Xbezier[i+1];
        Ybezier[i] -= gam[i] * Ybezier[i+1];
    }

    Xbezier[0] = Xspline[0];
    Ybezier[0] = Yspline[0];

    Xbezier[degree] = Xspline[degree];
    Ybezier[degree] = Yspline[degree];

    // draw the Splines
    for (i = 1; i <= degree; i++) {
        for (var j = 1; j < 5; j++) {
            var t  = 0.25 * j;
            var t1 = 1.0 - t;

            var xx = (t1 * t1 * t1) * (    Xspline[i-1]                 )
                   + (t1 * t1 * t ) * (2 * Xbezier[i-1] +     Xbezier[i])
                   + (t1 * t  * t ) * (    Xbezier[i-1] + 2 * Xbezier[i])
                   + (t  * t  * t ) * (                       Xspline[i]);
            var yy = (t1 * t1 * t1) * (    Yspline[i-1]                 )
                   + (t1 * t1 * t ) * (2 * Ybezier[i-1] +     Ybezier[i])
                   + (t1 * t  * t ) * (    Ybezier[i-1] + 2 * Ybezier[i])
                   + (t  * t  * t ) * (                       Yspline[i]);
            context.lineTo(xx,yy);
        }
    }
};


//
// get Point closest to the cursor
//
var sketcherGetClosestPoint = function () {
    // alert("sketcherGetClosestPoint()");

    var npnt  = sketch.pnt.x.length;
    var ibest = -1;
    var dbest = 999999;

    for (var ipnt = 0; ipnt < npnt; ipnt++) {
        var dtest = Math.abs(wv.cursorX - sketch.pnt.x[ipnt])
                  + Math.abs(wv.cursorY - sketch.pnt.y[ipnt]);
        if (dtest < dbest) {
            ibest = ipnt;
            dbest = dtest;
        }
    }

    return ibest;
};


//
// get Segment closest to the cursor
//
var sketcherGetClosestSegment = function () {
    // alert("sketcherGetClosestSegment()");

    var nseg  = sketch.seg.type.length;
    var ibest = -1;
    var dbest = 999999;

    for (var iseg = 0; iseg < nseg; iseg++) {
        var dtest = Math.abs(wv.cursorX - sketch.seg.xm[iseg])
                  + Math.abs(wv.cursorY - sketch.seg.ym[iseg]);
        if (dtest < dbest) {
            ibest = iseg;
            dbest = dtest;
        }
    }

    return ibest;
};


//
// get W or D constrint closest to the cursor
//
var sketcherGetClosestConstraint = function () {
    // alert("sketcherGetClosestConstraint()");

    var ncon  = sketch.con.type.length;
    var ibest = undefined;
    var dbest = 999999;

    for (var icon = 0; icon < ncon; icon++) {
        if        (sketch.con.type[icon] == "W") {
            var xmid = (  sketch.pnt.x[sketch.con.index1[icon]]
                        + sketch.pnt.x[sketch.con.index2[icon]]) / 2;
            var ymid = (2*sketch.pnt.y[sketch.con.index1[icon]]
                        + sketch.pnt.y[sketch.con.index2[icon]]) / 3;
        } else if (sketch.con.type[icon] == "D") {
            var xmid = (2*sketch.pnt.x[sketch.con.index1[icon]]
                        + sketch.pnt.x[sketch.con.index2[icon]]) / 3;
            var ymid = (  sketch.pnt.y[sketch.con.index1[icon]]
                        + sketch.pnt.y[sketch.con.index2[icon]]) / 2;
        } else {
            continue;
        }

        var dtest = Math.abs(wv.cursorX - xmid) + Math.abs(wv.cursorY - ymid);
        if (dtest < dbest) {
            ibest = icon;
            dbest = dtest;
        }
    }

    return ibest;
}
