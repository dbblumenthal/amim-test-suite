package de.mpg.mpiinf.ag1.kpm;

import java.awt.Component;
import java.awt.Font;
import java.awt.GridLayout;
import java.util.Vector;
import javax.swing.JFrame;
import javax.swing.JPanel;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.category.CategoryItemRenderer;
import org.jfree.chart.renderer.category.StatisticalBarRenderer;
import org.jfree.data.statistics.StatisticalCategoryDataset;
import org.jfree.data.xy.XYDataset;

/**
 * Window for creating and displaying graphs.
 *
 * @author timo
 */
public class MyWindow extends JFrame {

    /**
     * Creates a new Window for displaying components.
     *
     * @param title The frame title.
     * @param cs Components to display.
     */
    public MyWindow(String title, Vector<Component> cs) {
        super(title);
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        JPanel contentPane = new JPanel(new GridLayout(0, (int) Math.sqrt(cs.size())));

        for (Component c : cs) {
            contentPane.add(c);
        }
        setContentPane(contentPane);
        pack();
        setVisible(true);
    }

    /**
     * Creates a new Window for displaying components.
     *
     * @param title The frame title.
     * @param cs Components to display.
     * @param drawOnlyFirstComponent true: draw only the first component
     */
    public MyWindow(String title, Vector<Component> cs, boolean drawOnlyFirstComponent) {
        super(title);
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        JPanel contentPane = new JPanel(new GridLayout(0, (int) Math.sqrt(cs.size())));

        for (Component c : cs) {
            contentPane.add(c); break;
        }
        setContentPane(contentPane);
        pack();
        setVisible(true);
    }

//    /**
//     * Creates a new Window for displaying components.
//     *
//     * @param title The frame title.
//     * @param cs Components to display.
//     * @param axisStyle select scaling: 0 linear yAxis, 1 log yAxis, 2 log xAxis/yAxis
//     */
//    public MyWindow(String title, Vector<Component> cs, int axisStyle) {
//        super(title);
//        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
//
//        for (Component c : cs) {
//            if (axisStyle == 0) {
//                //do linear scale - here: nothing to change
//            } else if (axisStyle == 1) {
//                //do logarithmic scale
//                XYPlot xyplot = ((XYPlot) ((ChartPanel) c).getChart().getPlot());
//                String yAxis = xyplot.getRangeAxis().getLabel();
//                ValueAxis va = new LogarithmicAxis(yAxis);
////                ValueAxis va = new LogAxis(yAxis); //alternative zu LogarithmicAxis (looks different)
//                xyplot.setRangeAxis(va);
//            } else if (axisStyle == 2) {
//                //do logarithmic scale on xAxis and yAxis
//                XYPlot xyplot = ((XYPlot) ((ChartPanel) c).getChart().getPlot());
//                String yAxis = xyplot.getRangeAxis().getLabel();
//                String xAxis = xyplot.getDomainAxis().getLabel();
//                ValueAxis ya = new LogarithmicAxis(yAxis);
//                ValueAxis xa = new LogarithmicAxis(xAxis);
////                ValueAxis va = new LogAxis(yAxis); //alternative zu LogarithmicAxis (looks different)
//                xyplot.setRangeAxis(ya);
//                xyplot.setDomainAxis(xa);
//            }
//            break;
//            /* "important hack" - otherwise
//            java.lang.RuntimeException: Values less than or equal to zero not allowed with logarithmic axis
//             * (no effective fix found yet)
//             */
//        }
//
//        JPanel contentPane = new JPanel(new GridLayout(0, (int) Math.sqrt(cs.size())));
//
//        for (Component c : cs) {
//            contentPane.add(c);
//        }
//        setContentPane(contentPane);
//        pack();
//        setVisible(true);
//    }

    public static ChartPanel buildStatChart(StatisticalCategoryDataset dataset, String label) {
        return buildStatChart(dataset, label, "Number of Iterations");
    }

    public static ChartPanel buildStatChart(StatisticalCategoryDataset dataset, String label, String label2) {

        final CategoryAxis xAxis = new CategoryAxis("Algorithm");
        xAxis.setLowerMargin(0.01d); // percentage of space before first bar
        xAxis.setUpperMargin(0.01d); // percentage of space after last bar
        xAxis.setCategoryMargin(0.1d); // percentage of space between categories
        final ValueAxis yAxis = new NumberAxis(label2);
//        final ValueAxis yAxis = new LogarithmicAxis(label2);

        // define the plot
        final CategoryItemRenderer renderer = new StatisticalBarRenderer();
        final CategoryPlot plot = new CategoryPlot(dataset, xAxis, yAxis, renderer);

        final JFreeChart chart = new JFreeChart(label,
                new Font("Helvetica", Font.BOLD, 14),
                plot,
                true);
        /*        try {
        ChartUtilities.saveChartAsJPEG(new File("Out1.jpg"), chart, 500, 300);
        } catch (IOException ex) {
        Logger.getLogger(MyWindow.class.getName()).log(Level.SEVERE, null, ex);
        }
         */

        ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new java.awt.Dimension(500, 270));

        return chartPanel;
    }

    public static ChartPanel buildXYChart(XYDataset data, String caption, String xAxis, String yAxis) {
        JFreeChart chart = ChartFactory.createXYLineChart(caption, xAxis, yAxis,
                data, PlotOrientation.VERTICAL, true, true, false);

        ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new java.awt.Dimension(400, 300));

        return chartPanel;
    }

    public static ChartPanel buildXYChartLogarithmicY(XYDataset data, String caption, String xAxis, String yAxis) {
        JFreeChart chart = ChartFactory.createXYLineChart(caption, xAxis, yAxis,
                data, PlotOrientation.VERTICAL, true, true, false);
        ((XYPlot) (chart.getPlot())).setRangeAxis(new LogarithmicAxis(yAxis));

        ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new java.awt.Dimension(400, 300));

        return chartPanel;
    }
}
