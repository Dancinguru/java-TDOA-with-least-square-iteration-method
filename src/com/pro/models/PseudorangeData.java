package com.pro.models;

public class PseudorangeData {
    private double[] receiverPosition;
    private double pseudorange;
    private double error;

    public PseudorangeData(double[] receiverPosition, double pseudorange, double error) {
        this.receiverPosition = receiverPosition;
        this.pseudorange = pseudorange;
        this.error = error;
    }

    public double[] getReceiverPosition() {
        return receiverPosition;
    }

    public double getPseudorange() {
        return pseudorange;
    }

    public double getError() {
        return error;
    }
}
