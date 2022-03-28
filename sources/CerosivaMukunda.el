/*-----------------------------------------------------------------------------------------
 LIBRARY: LPRES_EXAMPLES_ALVAROG234
 FILE: CerosivaMukunda
 CREATION DATE: 28/03/2022
-----------------------------------------------------------------------------------------*/
FUNCTION REAL eta_mukunda(IN REAL g0,IN REAL gth,IN INTEGER nodos ,OUT REAL eta)
BODY

IF (eta<1) THEN
eta = 1
END IF

RETURN eta

END FUNCTION