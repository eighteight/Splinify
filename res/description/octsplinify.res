CONTAINER Octsplinify
{
	NAME Octsplinify;
	INCLUDE Obase;

	GROUP ID_OBJECTPROPERTIES
	{
		REAL CTTSPOBJECT_MAXSEG {MIN 0.; STEP 1.0;}
		BOOL CTTSPOBJECT_REL {}
	}

INCLUDE Osplineprimitive;
}
