package de.ovgu.mpa.validator;

import java.util.concurrent.atomic.AtomicLong;

public class ProgressTask implements Runnable {

    ProgressBar progressBar;
    AtomicLong done;
    long total;

    public ProgressTask(AtomicLong done, long total) {
        progressBar = new ProgressBar();
        this.done = done;
        this.total = total;
    }

    @Override
    public void run() {
        while (done.get() < total) {
            progressBar.update(done.get(), total);
            try {
                Thread.sleep(500);
            } catch (InterruptedException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }
        progressBar.update(done.get(), total);
    }
    
}