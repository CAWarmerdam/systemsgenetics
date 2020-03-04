package nl.systemsgenetics.pgsbasedmixupmapper.robustnessanalyzer;

public class RobustnessAnalyzerException extends Exception {
    public RobustnessAnalyzerException(String message, Exception e) {
        super(message, e);
    }
}
