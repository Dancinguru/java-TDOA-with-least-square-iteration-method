package com.pro;

import com.pro.Enum.Constants;
import com.pro.models.Measurement;
import com.pro.models.Receiver;
import com.pro.utils.Geodesy;
import com.pro.utils.Solver;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

public class Main {

    public static void main(String[] args) {
//        double[] rpos0 = Geodesy.llh2ecef(new double[]{51.45564, -0.96672, 99.9744});
//        double[] rpos1 = Geodesy.llh2ecef(new double[]{51.77872, -1.08092, 75.8952});
//        double[] rpos2 = Geodesy.llh2ecef(new double[]{51.42375, -1.09686, 48.1584});
//        double[] rpos3 = Geodesy.llh2ecef(new double[]{53.4288, -6.26331, 89.6112});

        double[] rpos0 = Geodesy.llh2ecef(new double[]{51.45564, -0.96672, 99.9744});
        double[] rpos1 = Geodesy.llh2ecef(new double[]{51.77872, -1.08092, 75.8952});
        double[] rpos2 = Geodesy.llh2ecef(new double[]{51.42375, -1.09686, 48.1584});
        double[] rpos3 = Geodesy.llh2ecef(new double[]{53.4288, -6.26331, 89.6112});
        double[] rpos4 = Geodesy.llh2ecef(new double[]{52.4288, -3.26331, 189.6112});
        double[] rpos5 = Geodesy.llh2ecef(new double[]{52.0288, -2.26331, 289.6112});
        double testalt = 9144;
        double timeError = 0.000001;

//        double[] testplane = Geodesy.llh2ecef(new double[]{51.50614956103109, -1.5365650867828802, testalt});
//        double t0 = 10 + Geodesy.norm(testplane, rpos0) / Constants.Cair + (0.5 - Math.random() * 1e-7);
//        double t1 = 10 + Geodesy.norm(testplane, rpos1) / Constants.Cair + (0.5 - Math.random() * 1e-7) ;
//        double t2 = 10 + Geodesy.norm(testplane, rpos2) / Constants.Cair + (0.5 - Math.random() * 1e-7) ;
//        double t3 = 10 + Geodesy.norm(testplane, rpos3) / Constants.Cair + (0.5 - Math.random() * 1e-7) ;

        double[] testplane = Geodesy.llh2ecef(new double[]{51.50614956103109, -1.5365650867828802, testalt});
        double t0 = 10 + Geodesy.norm(testplane, rpos0) / Constants.Cair + (0.5 - Math.random()) * 1e-7;
        double t1 = 10 + Geodesy.norm(testplane, rpos1) / Constants.Cair + (0.5 - Math.random()) * 1e-7;
        double t2 = 10 + Geodesy.norm(testplane, rpos2) / Constants.Cair + (0.5 - Math.random()) * 1e-7;
        double t3 = 10 + Geodesy.norm(testplane, rpos3) / Constants.Cair + (0.5 - Math.random()) * 1e-7;
        double t4 = 10 + Geodesy.norm(testplane, rpos4) / Constants.Cair + (0.5 - Math.random()) * 1e-7;
        double t5 = 10 + Geodesy.norm(testplane, rpos5) / Constants.Cair + (0.5 - Math.random()) * 1e-7;

        List<Measurement> measurements = new ArrayList<>();
        measurements.add(new Measurement(new Receiver(rpos0), t0, timeError));
        measurements.add(new Measurement(new Receiver(rpos1), t1, timeError));
        measurements.add(new Measurement(new Receiver(rpos2), t2, timeError));
        measurements.add(new Measurement(new Receiver(rpos3), t3, timeError));
        measurements.add(new Measurement(new Receiver(rpos4), t4, timeError));
        measurements.add(new Measurement(new Receiver(rpos5), t5, timeError));

        measurements.sort(Comparator.comparingDouble(Measurement::getTimestamp));
        Solver solver = new Solver();
        double[] result = solver.solve(measurements, testalt, 10.0, measurements.get(0).getReceiver().getPosition());

        if (result != null) {
            System.out.println("Position: " + Geodesy.ecef2llh(result)[0] + ", "+ Geodesy.ecef2llh(result)[1] + ", "+ Geodesy.ecef2llh(result)[2]);
            double error = Geodesy.norm(Geodesy.subtract(result, testplane));
            System.out.printf("Error: %.9fm\n", error);
        } else {
            System.out.println("Did not converge");
        }
    }
}