// ESP-plotter.js implements plotter functions for the Engineering Sketch Pad (ESP)
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

// interface functions that ESP-plotter.js provides
//    plotter.launch()
//    plotter.cmdUndo()
//    plotter.cmdSolve()
//    plotter.cmdSave()
//    plotter.cmdQuit()
//
//    plotter.cmdHome()                    not provided
//    plotter.cmdLeft()                    not provided
//    plotter.cmdRite()                    not provided
//    plotter.cmdBotm()                    not provided
//    plotter.cmdTop()                     not provided
//    plotter.cmdIn()                      not provided
//    plotter.cmdOut()                     not provided
//
//    plotter.mouseDown(e)                 not provided
//    plotter.mouseMove(e)                 not provided
//    plotter.mouseUp(e)                   not provided
//    plotter.mouseWheel(e)                not provided
//    plotter.mouseLeftCanvas(e)           not provided
//
//    plotter.keyPress(e)                  not provided
//    plotter.keyDown(e)                   not provided
//    plotter.keyUp(e)                     not provided
//    plotter.keyPressPart1(myKeyPress)    not provided
//    plotter.keyPressPart2(picking, gprim)not provided
//    plotter.updateKeyWindow()
//
//    plotter.timLoadCB(text)
//    plotter.timSaveCB(text)
//    plotter.timQuitCB(text)
//    plotter.timMesgCB(text)

"use strict";


//
// name of TIM object
//
plotter.name = "plotter";


//
// callback when Plotter is launched
//
plotter.launch = function () {
    alert("in plotter.launch()");

};


//
// callback when "Plotter->Undo" is pressed (called by ESP.html)
//
plotter.cmdUndo = function () {
    alert("plotter.cmdUndo() is not implemented");
};


//
// callback when "solveButton" is pressed
//
plotter.cmdSolve = function () {
    alert("in plotter.cmdSolve() is not implemented");

};


//
// callback when "Plotter->Save" is pressed (called by ESP.html)
//
plotter.cmdSave = function () {
    alert("in plotter.cmdSave() is not implemented");

};


//
// callback when "Plotter->Quit" is pressed (called by ESP.html)
//
plotter.cmdQuit = function () {
    alert("in plotter.cmdQuit() is not implemented");

};


//
// callback when "homeButton" is pressed (calles by ESP.html)
//
//plotter.cmdHome = function () {
//    main.cmdHome();
//};


//
// callback when "leftButton" is pressed (calles by ESP.html)
//
//plotter.cmdLeft = function () {
//    main.cmdLeft();
//};


//
// callback when "riteButton" is pressed (calles by ESP.html)
//
//plotter.cmdRite = function () {
//    main.cmdRite();
//};


//
// callback when "botmButton" is pressed (calles by ESP.html)
//
//plotter.cmdBotm = function () {
//    main.cmdBotm();
//};


//
// callback when "topButton" is pressed (calles by ESP.html)
//
//plotter.cmdTop = function () {
//    main.cmdTop();
//};


//
// callback when "inButton" is pressed (calles by ESP.html)
//
//plotter.cmdIn = function () {
//    main.cmdIn();
//};


//
// callback when "outButton" is pressed (calles by ESP.html)
//
//plotter.cmdOut = function () {
//    main.cmdOut();
//};


//
// callback when any mouse is pressed in canvas
//
//plotter.mouseDown = function(e) {
//    main.mouseDown(e);
//}


//
// callback when any mouse moves in canvas
//
//plotter.mouseMove = function(e) {
//    main.mouseMove(e);
//}


//
// callback when the mouse is released in canvas
//
//plotter.mouseUp = function(e) {
//    main.mouseUp(e);
//}


//
// callback when the mouse wheel is rolled in canvas
//
//plotter.mouseWheel = function(e) {
//    main.mouseWheel(e);
//}


//
// callback when the mouse leaves the canvas
//
//plotter.mouseLeftCanvas = function(e) {
//    main.mouseLeftCanvas(e);
//}


//
// callback when a key is pressed
//
//plotter.keyPress = function (e) {
//    main.keyPress(e);
//}


//
// callback when an arror... or shift is pressed (needed for Chrome)
//
//plotter.keyDown = function (e) {
//    main.keyDown(e);
//};


//
// callback when a shift is released (needed for Chrome)
//
//plotter.keyUp = function (e) {
//    main.keyUp(e);
//};


//
// callback for first part of a keypress that is not recognized by wvUpdateUI
//
//plotter.keyPressPart1 = function(myKeyPress) {
//    // alert("in plotter.keyPressPart1(myKeyPress="+myKeyPress+")");
//    return 0;
//};


//
// callback for second part of a keypress that is not recognized by wvUpdateUI
//
//plotter.keyPressPart2 = function(picking, gprim) {
//    // alert("in plotter.keyPressPart2(picking="+picking+"   gprim="+gprim+")");
//};


//
// function to update the key window
//
plotter.updateKeyWindow = function () {
    alert("in plotter.updateKeyWindow()");

};


//
// callback when timLoad returns
//
plotter.timLoadCB = function (text) {
    alert("in plotter.timLoadCB: "+text+" ignored)");

};


//
// callback when timSave returns
//
plotter.timSaveCB = function (text) {
    alert("in plotter.timSaveCB: "+text+" (ignored)");

};


//
// callback when timQuit returns
//
plotter.timQuitCB = function (text) {
    //alert("in plotter.timQuitCB: "+text+" (ignored)");

    var webViewer        = document.getElementById("WebViewer");
    var sketcherForm     = document.getElementById("sketcherForm");

    webView.hidden      = false;
    sketcherForm.hidden = true;
};


//
// callback when timMesg returns
//
plotter.timMesgCB = function (text) {
    //alert("in plotter.timMesgCB: "+text);

    if (text.substring(0, 13) != "plotter|show|") {
        return;
    }

    // unwrap the JSON data
    var data = JSON.parse(text.substring(13,text.length-1));

    // show the sketcher (canvas)
    var webViewer        = document.getElementById("WebViewer");
    var sketcherForm     = document.getElementById("sketcherForm");

    webViewer.hidden    = true;
    sketcherForm.hidden = false;

    // set up graphics context
    var canvas  = document.getElementById("sketcher");
    var context = canvas.getContext("2d");

    // clear the screen
    context.clearRect(0, 0, canvas.width, canvas.height);

    context.lineCap  = "round";
    context.lineJoin = "round";

    // find the extrema of the data
    var xmin = data.lines[0].x[0];
    var xmax = data.lines[0].x[0];
    var ymin = data.lines[0].y[0];
    var ymax = data.lines[0].y[0];

    for (var iline = 0; iline < data.lines.length; iline++) {
        for (var ipnt = 0; ipnt < data.lines[iline].x.length; ipnt++) {
            if (xmin > data.lines[iline].x[ipnt]) {
                xmin = data.lines[iline].x[ipnt];
            }
            if (xmax < data.lines[iline].x[ipnt]) {
                xmax = data.lines[iline].x[ipnt];
            }
            if (ymin > data.lines[iline].y[ipnt]) {
                ymin = data.lines[iline].y[ipnt];
            }
            if (ymax < data.lines[iline].y[ipnt]) {
                ymax = data.lines[iline].y[ipnt];
            }
        }
    }

    if        (Math.abs(xmax-xmin) < 1e-6) {
        alert("all x="+(xmin+xmax)/2);
        return;
    } else if (Math.abs(ymax-ymin) < 1e-6) {
        alert("all y="+(ymin+ymax)/2);
        return;
    }

    // compute the scale factors
    var ixmin =                100;
    var ixmax = canvas.width  - 20;
    var iymin = canvas.height - 50;
    var iymax =                 20;

    var xm = (ixmax - ixmin) / (xmax - xmin);
    var xa =  ixmin - xmin * xm;
    var ym = (iymin - iymax) / (ymin - ymax);
    var ya =  iymax - ymax * ym;

    // draw axes
    context.lineWidth = 1;
    context.beginPath();
    context.moveTo(ixmin, iymax);
    context.lineTo(ixmin, iymin);
    context.lineTo(ixmax, iymin);
    context.stroke();

    context.font         = "12px Verdana";
    context.textBaseline = "middle";
    context.fillStyle    = "black";

    // title, xlabel, and x-axis values
    context.textAlign    = "center";
    context.fillText(data.title,  (ixmin+ixmax)/2, iymax-10);
    context.fillText(data.xlabel, (ixmin+ixmax)/2, iymin+35);

    for (var ii = 0; ii < 5; ii++) {
        var ix = ixmin + ii/4 * (ixmax - ixmin);
        var xx =  xmin + ii/4 * ( xmax -  xmin);
        context.fillText(xx.toPrecision(3), ix, iymin+20);
    }

    // ylabel
    context.textAlign = "left";
    context.fillText(data.ylabel, 10, (iymin+iymax)/2);

    // y-axis values
    context.textAlign    = "right";

    for (var jj = 0; jj < 4; jj++) {
        var iy = iymin + jj/3 * (iymax - iymin);
        var yy =  ymin + jj/3 * ( ymax -  ymin);
        context.fillText(yy.toPrecision(3), ixmin-10, iy);
    }

    // draw the lines
    context.lineWidth = 3;
    for (iline = 0; iline < data.lines.length; iline++) {

        // get color
        if        (data.lines[iline].style.indexOf("r") >= 0) {     // red
            context.strokeStyle = "red";
            context.fillStyle   = "red";
        } else if (data.lines[iline].style.indexOf("g") >= 0) {     // green
            context.strokeStyle = "green";
            context.fillStyle   = "green";
        } else if (data.lines[iline].style.indexOf("b") >= 0) {     // blue
            context.strokeStyle = "blue";
            context.fillStyle   = "blue";
        } else if (data.lines[iline].style.indexOf("c") >= 0) {     // cyan
            context.strokeStyle = "cyan";
            context.fillStyle   = "cyan";
        } else if (data.lines[iline].style.indexOf("m") >= 0) {     // magenta
            context.strokeStyle = "magenta";
            context.fillStyle   = "magenta";
        } else if (data.lines[iline].style.indexOf("y") >= 0) {     // yellow
            context.strokeStyle = "yellow";
            context.fillStyle   = "yellow";
        } else if (data.lines[iline].style.indexOf("k") >= 0) {     // black
            context.strokeStyle = "black";
            context.fillStyle   = "black";
        } else if (data.lines[iline].style.indexOf("w") >= 0) {     // white
            context.strokeStyle = "white";
            context.fillStyle   = "white";
        }

        // get line styles
        var showLine = 0;
        if        (data.lines[iline].style.indexOf("-") >= 0) {     // solid
            context.setLineDash([]);
            showLine++;
        } else if (data.lines[iline].style.indexOf("_") >= 0) {     // dashed
            context.setLineDash([20, 5]);
            showLine++;
        } else if (data.lines[iline].style.indexOf(":") >= 0) {     // dotted
            context.setLineDash([2,5]);
            showLine++;
        } else if (data.lines[iline].style.indexOf(";") >= 0) {     // dot-dash
            context.setLineDash([20,5, 2, 5]);
            showLine++;
        }

        // draw the lines
        if (showLine > 0) {
            context.beginPath();
            var x = xa + xm * data.lines[iline].x[0];
            var y = ya + ym * data.lines[iline].y[0];
            context.moveTo(x, y);

            for (var ipnt = 1; ipnt < data.lines[iline].x.length; ipnt++) {
                x = xa + xm * data.lines[iline].x[ipnt];
                y = ya + ym * data.lines[iline].y[ipnt];
                context.lineTo(x, y);
            }
            context.stroke();
        }
        context.setLineDash([]);

        // add the symbols
        if        (data.lines[iline].style.indexOf("o") >= 0) {     // circle
            for (ipnt = 0; ipnt < data.lines[iline].x.length; ipnt++) {
                x = xa + xm * data.lines[iline].x[ipnt];
                y = ya + ym * data.lines[iline].y[ipnt];
                context.beginPath();
                context.arc(x, y, 5, 0, 2*Math.PI);
                context.fill();
            }
        } else if (data.lines[iline].style.indexOf("x") >= 0) {     // X
            for (ipnt = 0; ipnt < data.lines[iline].x.length; ipnt++) {
                x = xa + xm * data.lines[iline].x[ipnt];
                y = ya + ym * data.lines[iline].y[ipnt];
                context.beginPath();
                context.moveTo(x-5, y-5);
                context.lineTo(x+5, y+5);
                context.stroke();

                context.beginPath();
                context.moveTo(x-5, y+5);
                context.lineTo(x+5, y-5);
                context.stroke();
            }
        } else if (data.lines[iline].style.indexOf("+") >= 0) {     // +
            for (ipnt = 0; ipnt < data.lines[iline].x.length; ipnt++) {
                x = xa + xm * data.lines[iline].x[ipnt];
                y = ya + ym * data.lines[iline].y[ipnt];
                context.beginPath();
                context.moveTo(x, y-5);
                context.lineTo(x, y+5);
                context.stroke();

                context.beginPath();
                context.moveTo(x-5, y);
                context.lineTo(x+5, y);
                context.stroke();
            }
        } else if (data.lines[iline].style.indexOf("s") >= 0) {     // square
            for (ipnt = 0; ipnt < data.lines[iline].x.length; ipnt++) {
                x = xa + xm * data.lines[iline].x[ipnt];
                y = ya + ym * data.lines[iline].y[ipnt];
                context.fillRect(x-5, y-5, 10, 10);
            }
        }
    }

    context.strokeStyle = "black";
    context.fillStyle   = "black";
};


// /////////////////////////////////////////////////////////////////////

