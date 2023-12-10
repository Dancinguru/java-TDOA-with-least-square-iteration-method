package com.pro.models;

public class Measurement {
    private Receiver receiver;
    private double timestamp;
    private double variance;

    public Measurement(Receiver receiver, double timestamp, double variance) {
        this.receiver = receiver;
        this.timestamp = timestamp;
        this.variance = variance;
    }

    public Receiver getReceiver() {
        return receiver;
    }

    public double getTimestamp() {
        return timestamp;
    }

    public double getVariance() {
        return variance;
    }
}
