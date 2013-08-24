CONTAINER Octsplinify
{
	NAME Octsplinify;
	INCLUDE Obase;

	GROUP ID_OBJECTPROPERTIES
	{
		REAL CTTSPOBJECT_MAXSEG { UNIT METER;	MIN 0.0; CUSTOMGUI REALSLIDER; STEP 0.01; MINSLIDER 0.0; MAXSLIDER 100.0;}
		LONG START_FRAME {MIN 0; CUSTOMGUI LONGSLIDER; STEP 1; MINSLIDER 0; MAXSLIDER 10000;}
		LONG END_FRAME {MIN 0; CUSTOMGUI LONGSLIDER; STEP 1; MINSLIDER 0; MAXSLIDER 10000;}
		LINK CTT_OBJECT_LINK { ACCEPT { Obase; } }
		LONG CTTSPOBJECT_WINDOW { MIN 1; MAX 1000; CUSTOMGUI LONGSLIDER; }
		LONG CTT_SPLINE_PERCENTAGE { MIN 1; MAX 100; CUSTOMGUI LONGSLIDER; }
		BOOL CTTSPOBJECT_REL {}
	}

INCLUDE Osplineprimitive;
}
