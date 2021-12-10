			CAPS Prerequisite -- UDUNITS2

	CAPS requires the use of the Open Source Project UDUNITS2 for all
unit conversions. Since there are no prebuilt package distributions for
the Mac and Windows, the CAPS build procedure copies the DLL/DyLibs from
this directory to the lib directory of ESP. Because most LINUX distributions
contain a UDUNITS2 package, this is not done under LINUX and that package
must be installed by the system procedure (based on the LINUX flavor). This
is done to avoid having two differing versions on the same system (which
would only cause problems).
