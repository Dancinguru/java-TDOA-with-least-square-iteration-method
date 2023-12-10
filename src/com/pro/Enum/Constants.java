package com.pro.Enum;

public  class Constants {
    public static double Cair = 299792458 / 1.0003;

    public static double DTOR = Math.PI / 180.0;
    public static double RTOD = 180.0 / Math.PI;

    public static double FTOM = 0.3038;
    public static double MTOF = 1.0/FTOM;

    public static double MS_TO_KTS = 1.9438;

    public static double MS_TO_FPM = MTOF * 60;
    public static final double MAX_RANGE = 1.0e6; // Adjust the maximum range as needed
    public static final int SOLVER_MAXFEV = 100; // Adjust the maximum number of function evaluations as needed

    public static double[][] teststations1 = {{51.45564, -0.96672, 99.9744}, {51.77872, -1.08092, 75.8952},
            {51.42375, -1.09686, 48.1584}, {53.4288, -6.26331, 89.6112}};
    public static int testalt1 = 9144;
    public static double[] testplane1 = new double[]{51.50614956103109, -1.5365650867828802, testalt1};
    public static double[][] teststations2 = {{51.45564, -0.96672, 99.9744}, {51.77872, -1.08092, 75.8952},
            {51.42375, -1.09686, 48.1584}, {53.4288, -6.26331, 89.6112}};
    public static double[] testplane2 = new double[]{50.50614956103109, 1.5365650867828802, testalt1};
    public static double[][] teststations3 =  {{51.45564, -0.96672, 99.9744}, {51.77872, -1.08092, 75.8952},
            {51.92375, -1.09686, 48.1584}};
    public static double[] testplane3 = new double[]{51.50614956103109, -1.0365650867828802, testalt1};
}
