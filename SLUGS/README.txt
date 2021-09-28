Tutorial1 test case
===================

   serveCSM ../data/Slugs/tutorial1_setup -batch -egg PDT2
   Slugs    ../data/Slugs/tutorial1_setup -batch -jrnl ../data/Slugs/tutorial1_setup.jrnl
   Slugs    ../data/Slugs/tutorial1       -batch -jrnl ../data/Slugs/tutorial1.jrnl
   serveCSM ../data/Slugs/tutorial1       -batch

Cylinder test case
==================

   serveCSM ../data/Slugs/cylinder_0 -batch
   Slugs    ../data/Slugs/cylinder_0 -batch -jrnl ../data/Slugs/cylinder_0.jrnl
   serveCSM ../data/Slugs/cylinder_1

Blend test case
================

   serveCSM ../data/Slugs/blend12_0 -batch
   Slugs    ../data/Slugs/blend12_0 -batch -jrnl ../data/Slugs/blend12_0.jrnl
   serveCSM ../data/Slugs/blend12_1

Wing/body test case
===================

   serveCSM ../data/Slugs/wingBody_0 -batch
   Slugs    ../data/Slugs/wingBody_0 -batch -jrnl ../data/Slugs/wingBody_0.jrnl
   Slugs    ../data/Slugs/wingBody_1 -batch -jrnl ../data/Slugs/wingBody_1.jrnl
   serveCSM ../data/Slugs/wingBody_2

Demo1 test case
===============

   serveCSM ../data/Slugs/demo1_0 -batch -egg PDT2
   Slugs    ../data/Slugs/demo1_0 -batch -jrnl ../data/Slugs/demo1_0.jrnl
   serveCSM ../data/Slugs/demo1_1

Damaged test case
=================

   serveCSM ../data/Slugs/damaged_0 -batch
   Slugs    ../data/Slugs/damaged_0 -batch -jrnl ../data/Slugs/damaged_0.jrnl
   Slugs    ../data/Slugs/damaged_1 -batch -jrnl ../data/Slugs/damaged_1.jrnl
   serveCSM ../data/Slugs/damaged_2

myGlider test case
==================

   yet to come

myPlane test case
=================

   yet to come


