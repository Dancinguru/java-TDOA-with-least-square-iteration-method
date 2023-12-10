package com.pro.utils;

import com.pro.Enum.Constants;
import com.pro.models.Measurement;
import com.pro.models.PseudorangeData;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import static com.pro.Enum.Constants.MAX_RANGE;
import static com.pro.Enum.Constants.SOLVER_MAXFEV;

public class Geodesy {
    public static final double WGS84_A = 6378137.0;
    private static final double WGS84_F = 1.0 / 298.257223563;
    private static final double WGS84_B = WGS84_A * (1 - WGS84_F);
    public static final double WGS84_ECC_SQ = 1 - WGS84_B * WGS84_B / (WGS84_A * WGS84_A);
    private static final double WGS84_ECC = Math.sqrt(WGS84_ECC_SQ);

    // Average radius for a spherical Earth
    private static final double SPHERICAL_R = 6371e3;

    // Some derived values
    private static final double _wgs84_ep = Math.sqrt((WGS84_A * WGS84_A - WGS84_B * WGS84_B) / (WGS84_B * WGS84_B));
    private static final double _wgs84_ep2_b = _wgs84_ep * _wgs84_ep * WGS84_B;
    private static final double _wgs84_e2_a = WGS84_ECC_SQ * WGS84_A;

    public static double[] llh2ecef(double[] llh) {
        double lat = llh[0] * Constants.DTOR;
        double lng = llh[1] * Constants.DTOR;
        double alt = llh[2];

        double slat = Math.sin(lat);
        double slng = Math.sin(lng);
        double clat = Math.cos(lat);
        double clng = Math.cos(lng);

        double d = Math.sqrt(1 - (slat * slat * WGS84_ECC_SQ));
        double rn = WGS84_A / d;

        double x = (rn + alt) * clat * clng;
        double y = (rn + alt) * clat * slng;
        double z = (rn * (1 - WGS84_ECC_SQ) + alt) * slat;

        return new double[]{x, y, z};
    }

    public static double[] ecef2llh(double[] ecef) {
        double x = ecef[0];
        double y = ecef[1];
        double z = ecef[2];

        double lon = Math.atan2(y, x);

        double p = Math.sqrt(x * x + y * y);

        double th = Math.atan2(WGS84_A * z, WGS84_B * p);

        double lat = Math.atan2(z + _wgs84_ep2_b * Math.sin(th) * Math.sin(th) * Math.sin(th),
                p - _wgs84_e2_a * Math.cos(th) * Math.cos(th) * Math.cos(th));

        double N = WGS84_A / Math.sqrt(1 - WGS84_ECC_SQ * Math.sin(lat) * Math.sin(lat));
        double alt = p / Math.cos(lat) - N;

        return new double[]{lat * Constants.RTOD, lon * Constants.RTOD, alt};
    }
    public static double greatCircle(double[] p0, double[] p1) {
        // Returns a great-circle distance in metres between two LLH points,
        // _assuming spherical earth_ and _ignoring altitude_. Don't use this if you
        // need a distance accurate to better than 1%.
        double lat0 = p0[0] * Constants.DTOR;
        double lon0 = p0[1] * Constants.DTOR;
        double lat1 = p1[0] * Constants.DTOR;
        double lon1 = p1[1] * Constants.DTOR;
        return SPHERICAL_R * Math.acos(
                Math.sin(lat0) * Math.sin(lat1) +
                        Math.cos(lat0) * Math.cos(lat1) * Math.cos(Math.abs(lon0 - lon1)));
    }

    public static double ecefDistance(double[] p0, double[] p1) {
        // Returns the straight-line distance in metres between two ECEF points.
        int len = Math.min(p0.length, p1.length);

        double sum = 0;
        for(int i = 0; i < len; i ++ ){
            double temp = p0[i] - p1[i];
            sum += temp * temp;
        }

        return Math.sqrt(sum);
    }

    public static double norm( final double[] vector )
    {
        double sum = 0;
        for ( final double value : vector )
        {
            sum += value * value;
        }
        return Math.sqrt( sum );
    }

    public static double norm(double[] v1, double[] v2) {
        int len = Math.min(v1.length, v2.length);

        double sum = 0;
        for(int i = 0; i < len; i ++ ){
            double temp = v1[i] - v2[i];
            sum += temp * temp;
        }
        return Math.sqrt(sum);
    }

    public static double[] subtract( final double[] vector1, final double[] vector2 )
    {
        final double[] result = new double[Math.min(vector1.length, vector2.length)];
        for ( int i = 0; i < Math.min(vector1.length, vector2.length); i++ )
        {
            result[i] = vector1[i] - vector2[i];
        }
        return result;
    }
}
