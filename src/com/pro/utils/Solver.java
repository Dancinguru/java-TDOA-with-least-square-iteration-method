package com.pro.utils;

import com.pro.Enum.Constants;
import com.pro.models.Measurement;
import com.pro.models.PseudorangeData;
import com.pro.models.Receiver;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.fitting.leastsquares.*;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleVectorValueChecker;
import org.apache.commons.math3.util.Pair;

import java.sql.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;

import static com.pro.utils.Geodesy.WGS84_A;
import static com.pro.utils.Geodesy.WGS84_ECC_SQ;

public class Solver {
    private static final Logger logger = Logger.getLogger(Solver.class.getName());

    private static final double MAX_RANGE = 500e3;
    private static final int SOLVER_MAXFEV = 50;
    private List<PseudorangeData> pseudorangeData = new ArrayList<>();

    private List<Receiver> receivers = new ArrayList<>();

    MultivariateJacobianFunction distancesToCurrentCenter = new MultivariateJacobianFunction() {

        public Pair<RealVector, RealMatrix> value(final RealVector point) {
            int len = pseudorangeData.size() + 1;
            RealVector value = new ArrayRealVector(len);
            RealMatrix jacobian = new Array2DRowRealMatrix(len, 4);
            double x = point.getEntry(0);
            double y = point.getEntry(1);
            double z = point.getEntry(2);

            for (int i = 0; i < pseudorangeData.size(); ++i) {
                PseudorangeData p = pseudorangeData.get(i);
                double[] pcenter = new double[]{point.getEntry(0), point.getEntry(1), point.getEntry(2)};
                double modelI = (Geodesy.ecefDistance(p.getReceiverPosition(), pcenter) - point.getEntry(3));
                value.setEntry(i, modelI);
                // Calculate the Jacobian entries
                double dx = p.getReceiverPosition()[0] - pcenter[0];
                double dy = p.getReceiverPosition()[1] - pcenter[1];
                double dz = p.getReceiverPosition()[2] - pcenter[2];
                double divisor = Geodesy.ecefDistance(p.getReceiverPosition(), pcenter);
                if(modelI == 0){
                    jacobian.setEntry(i, 0, 0);
                    jacobian.setEntry(i, 1, 0);
                    jacobian.setEntry(i, 2, 0);
                    jacobian.setEntry(i, 3, -1);
                }
                else {
                    jacobian.setEntry(i, 0, -dx / divisor);
                    jacobian.setEntry(i, 1, -dy / divisor);
                    jacobian.setEntry(i, 2, -dz / divisor);
                    jacobian.setEntry(i, 3, -1);
                }
            }
            double[] llh = Geodesy.ecef2llh(new double[]{point.getEntry(0), point.getEntry(1), point.getEntry(2)});
            double altitude_guess = llh[2];
            double[] altJacob = Geodesy.getJacobOfAlt(new double[]{point.getEntry(0), point.getEntry(1), point.getEntry(2)});
            double alt_guess = altitude_guess ;
            value.setEntry(pseudorangeData.size(), alt_guess);

            jacobian.setEntry(pseudorangeData.size(), 0, altJacob[0]);
            jacobian.setEntry(pseudorangeData.size(), 1, altJacob[1]);
            jacobian.setEntry(pseudorangeData.size(), 2, altJacob[2]);
            jacobian.setEntry(pseudorangeData.size(), 3, 0);

            return new Pair<RealVector, RealMatrix>(value,  jacobian);
        }
    };
    public double[] solve(List<Measurement> measurements, Double altitude, Double altitudeError, double[] initialGuess) {
        if (measurements.size() + (altitude == null ? 0 : 1) < 4) {
            throw new IllegalArgumentException("Not enough measurements available");
        }
        List<Double> distances = new ArrayList<Double>();

        double baseTimestamp = measurements.get(0).getTimestamp();
        for (Measurement measurement : measurements) {
            double pseudorange = (measurement.getTimestamp() - baseTimestamp) * Constants.Cair;
            double error = Math.sqrt(measurement.getVariance()) * Constants.Cair;
            pseudorangeData.add(new PseudorangeData(measurement.getReceiver().getPosition(), pseudorange, error));
            distances.add(pseudorange);
        }
        distances.add(altitude);
        double[] xGuess = new double[] { initialGuess[0], initialGuess[1], initialGuess[2], 0.0 };

        double[] prescribedDistances = new double[distances.size()];


        for(int i = 0; i < distances.size(); i ++ ){
            prescribedDistances[i] = distances.get(i);
        }

        LeastSquaresProblem problem = new LeastSquaresBuilder().
                start(xGuess).
                model(distancesToCurrentCenter).
                target(prescribedDistances).
                lazyEvaluation(false).
                maxEvaluations(50).
                maxIterations(1000).build();

        LeastSquaresOptimizer.Optimum optimum = new LevenbergMarquardtOptimizer().optimize(problem);

        double[] xEst = optimum.getPoint().toArray();
        return xEst;
    }

}
