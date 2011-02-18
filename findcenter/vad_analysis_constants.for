C  Physical constants.
      REAL EARTH_RADIUS,DLNRIDZ,PI,RAD_CONVERSION
      PARAMETER(EARTH_RADIUS=6371229.,
     $DLNRIDZ=-4.E-8,
     $PI=3.141592654,
     $RAD_CONVERSION=PI/180.)

C  Logical input and output unit numbers
      INTEGER USER_INPUT_UNIT,DATA_INPUT_UNIT,BINARY_OUTPUT_UNIT,
     $RESULTS_OUTPUT_UNIT,PLOT_OUTPUT_UNIT,LOG_OUTPUT_UNIT,
     $USER_INFO_UNIT
      PARAMETER(USER_INPUT_UNIT=1,
     $RESULTS_OUTPUT_UNIT=6,
     $DATA_INPUT_UNIT=8,
     $BINARY_OUTPUT_UNIT=7,
     $LOG_OUTPUT_UNIT=9,
     $PLOT_OUTPUT_UNIT=20,
     $USER_INFO_UNIT=21)

C  ELEVATION_STOP_CODE is written in the binary output file to
C  mark the end of a constant-elevation sweep.
C  VOLUME_STOP_CODE is written in the binary output file to mark the end
C  of the volume.
      INTEGER ELEVATION_STOP_CODE,VOLUME_STOP_CODE
      PARAMETER(ELEVATION_STOP_CODE=-9998,
     $VOLUME_STOP_CODE=-9999)

C  BADFLAG represents missing data in the program.
      REAL BADFLAG
      PARAMETER(BADFLAG=-32768.)

C  As gates are processed, they are assigned status codes.
C  Zero is the code for gates with good data.
C  Positive codes designate gates that are eliminated from further
C  calculations.
C  Negative codes designate gates whose data have been modified but
C  that are still used in the analysis.  
C  DEALIAS_CODE is the code for gates whose velocities have been dealiased.
C  GOODDATA_CODE is the code for gates with good and unmodified data.
C  BADELEV_CODE is the code for gates from rays with deviant elevations.
C  LOWVELO_CODE is the code for gates that have low velocity and might be
C  ground clutter.
C  LOWPWR_CODE is the code for gates that have low power.
C  LOWDZ_CODE is the code for gates that have low reflectivity factor.
C  LOWNC_CODE is the code for gates that have low ncp value.
C  CLEAN_CODE is the code for gates that were thrown out because their
C  velocity deviated from the fit by more than a specified number of
C  standard deviations.
C  RESID_CODE is the code for gates that were thrown out because their
C  velocity deviated from the fit by more than a specified amount.
C  BADFLG_CODE is the code for gates with no data.
      INTEGER DEALIAS_CODE,LOWVELO_CODE,LOWPWR_CODE,LOWDZ_CODE,
     $CLEAN_CODE,RESID_CODE,BADFLG_CODE,GOODDATA_CODE,BADELEV_CODE,
     $LOWNC_CODE
      PARAMETER(DEALIAS_CODE=-1,
     $GOODDATA_CODE=0,
     $BADELEV_CODE=1,
     $LOWVELO_CODE=2,
     $LOWPWR_CODE=3,
     $LOWDZ_CODE=4,
     $CLEAN_CODE=5,
     $RESID_CODE=6,
     $BADFLG_CODE=7,
     $LOWNC_CODE=8)

C  Plot parameters.
C  H1_1, H2_1, H1_2, H2_2, H1_3, H2_3, H1_4, H2_4, V1_1, V2_1, V1_2,
C  and V2_2 are the positions of the edges of the plotting rectangles
C  as fractions of the total plotting region.
C  LINE_SIZE controls the line width.
C  DASH_CODE_SOLID controls the style of solid lines.
C  SMALL_TIC and LARGE_TIC are the sizes of small and large axis tic marks
C  as fractions of the total plotting region .
C  NUM_OFFSET_LR is the distance between the left or right axis and the
C  edge of the axis numbers as a fraction of the horizontal length of
C  the plotting rectangle.
C  NUM_OFFSET_BT is the distance between the bottom or top axis and the
C  center of the axis numbers as a fraction of the vertical length of
C  the plotting rectangle.
C  TITLE_POSITION1 is the position of the center of the lowest line
C  of the plot title as a fraction of the vertical length of the total
C  plotting region.
C  TITLE_POSITION2 is the position of the center of the second lowest line
C  of the plot title as a fraction of the vertical length of the total
C  plotting region.
C  PLOT_RES is the exponent of 2 that specifies the plot resolution (the
C  number of plot increments across the total plotting region).
C  PLOT_RES_REF is the exponent of 2 that specifies a reference plot
C  resolution.  Character widths are specified in terms of plot
C  increments for this reference resolution and then converted
C  to the actual resolution used.
C  NUM_SIZE is the width of characters in axis numbers.
C  LABEL_SIZE is the width of characters in axis labels.
C  TITLE_SIZE is the width of characters in plot titles.
C  SYMBOL_SIZE is the width of characters in plotted symbols.
      REAL H1_1,H2_1,H1_2,H2_2,H1_3,H2_3,H1_4,H2_4,V1_1,V2_1,V1_2,V2_2
      PARAMETER(H1_1=0.0342,
     $H2_1=0.2402,
     $H1_2=0.2822,
     $H2_2=0.4882,
     $H1_3=0.5303,
     $H2_3=0.7363,
     $H1_4=0.7783,
     $H2_4=0.9843,
     $V1_1=0.4883,
     $V2_1=0.8877,
     $V1_2=0.0225,
     $V2_2=0.4219)
      INTEGER LINE_SIZE
      PARAMETER(LINE_SIZE=8)
      INTEGER DASH_CODE_SOLID
      PARAMETER(DASH_CODE_SOLID=65535)
      REAL SMALL_TIC,LARGE_TIC
      PARAMETER(SMALL_TIC=0.006,
     $LARGE_TIC=0.01)
      REAL NUM_OFFSET_LR,NUM_OFFSET_BT,TITLE_POSITION1,
     $TITLE_POSITION2,PANEL_TITLE_OFFSET
      PARAMETER(NUM_OFFSET_LR=0.01,
     $NUM_OFFSET_BT=0.03,
     $TITLE_POSITION1=0.9365,
     $TITLE_POSITION2=0.9658,
     $PANEL_TITLE_OFFSET=0.03)
      INTEGER PLOT_RES,PLOT_RES_REF
      PARAMETER(PLOT_RES=10,
     $PLOT_RES_REF=10)
      INTEGER NUM_SIZE,TITLE_SIZE,SYMBOL_SIZE,PANEL_TITLE_SIZE
      PARAMETER(NUM_SIZE=(8*2**PLOT_RES)/2**PLOT_RES_REF,
     $TITLE_SIZE=(16*2**PLOT_RES)/2**PLOT_RES_REF,
     $SYMBOL_SIZE=(6*2**PLOT_RES)/2**PLOT_RES_REF,
     $PANEL_TITLE_SIZE=(8*2**PLOT_RES)/2**PLOT_RES_REF)
