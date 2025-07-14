// templates/spectrum_lookup.js

// Add this debounce function at the top of your file to limit update frequency
function debounce(func, wait) {
    let timeout;
    return function(...args) {
        clearTimeout(timeout);
        timeout = setTimeout(() => func.apply(this, args), wait);
    };
}

// Function to find closest compound spectrum to a given retention time
function findClosestSpectrum(rt) {
    if (!massSpectra || massSpectra.length === 0) return null;
    
    let closestSpectrum = null;
    let minDistance = Infinity;
    
    for (let i = 0; i < massSpectra.length; i++) {
        const spectrum = massSpectra[i];
        const distance = Math.abs(spectrum.rt - rt);
        
        if (distance < minDistance) {
            minDistance = distance;
            closestSpectrum = spectrum;
        }
    }
    
    return closestSpectrum;
}

// Create a spectrum plot
function createSpectrumPlot(spectrumData) {
    const mzDiv = document.getElementById('mz-spectrum-display');
    mzDiv.innerHTML = '';
    
    const title = document.createElement('h3');
    title.textContent = `Mass Spectrum: ${spectrumData.compound_type.charAt(0).toUpperCase() + spectrumData.compound_type.slice(1)}: ${spectrumData.name} (RT: ${spectrumData.rt.toFixed(2)})`;
    title.style.color = spectrumData.color;
    mzDiv.appendChild(title);
    
    // Create a basic plot using Plotly
    const plotDiv = document.createElement('div');
    plotDiv.style.width = '100%';
    plotDiv.style.height = '300px';
    mzDiv.appendChild(plotDiv);
    
    // FIXED LOLLIPOP IMPLEMENTATION - using null values to break lines between stems
    var stemX = [];
    var stemY = [];
    
    for (let i = 0; i < spectrumData.mz.length; i++) {
        // Base of the stem
        stemX.push(spectrumData.mz[i]);
        stemY.push(0);
        
        // Top of the stem
        stemX.push(spectrumData.mz[i]);
        stemY.push(spectrumData.intensity[i]);
        
        // Add null values to break the line between stems (creates individual lollipops)
        if (i < spectrumData.mz.length - 1) {
            stemX.push(null);
            stemY.push(null);
        }
    }
    
    // Stem lines
    const stemTrace = {
        x: stemX,
        y: stemY,
        mode: 'lines',
        line: {
            color: spectrumData.color,
            width: 1.5
        },
        hoverinfo: 'none',
        showlegend: false
    };
    
    // Points at stem tops
    const pointsTrace = {
        x: spectrumData.mz,
        y: spectrumData.intensity,
        mode: 'markers',  // Removed text to avoid clutter
        marker: {
            color: spectrumData.color,
            size: 5
        },
        hovertemplate: 'm/z: %{x:.2f}<br>Intensity: %{y}<extra></extra>',
        name: 'm/z peaks',
        showlegend: false
    };
    
    Plotly.newPlot(plotDiv, [stemTrace, pointsTrace], {
        title: {
            text: 'Mass Spectrum',
            font: {
                size: 16
            }
        },
        xaxis: {
            title: {
                text: 'm/z',
                font: {
                    size: 14
                }
            }
        },
        yaxis: {
            title: {
                text: 'Relative Intensity',
                font: {
                    size: 14
                }
            }
        },
        margin: {
            t: 30,
            b: 50,
            l: 60,
            r: 10
        }
    });
}

// Set up event handlers after page load
document.addEventListener('DOMContentLoaded', function() {
    var myPlot = document.getElementById('main-plot');
    var lastSpectrum = null;
    var crosshairLine = null;
    
    // Create crosshair line
    function updateCrosshair(x) {
        if (crosshairLine) {
            // Update existing line
            Plotly.relayout(myPlot, {
                'shapes[0].x0': x,
                'shapes[0].x1': x
            });
        } else {
            // Create new line on first hover
            const update = {
                shapes: [{
                    type: 'line',
                    x0: x,
                    y0: 0,
                    x1: x,
                    y1: 1,
                    yref: 'paper',
                    line: {
                        color: 'rgba(255, 0, 0, 0.7)',
                        width: 1.5,
                        dash: 'solid'
                    }
                }]
            };
            Plotly.relayout(myPlot, update);
            crosshairLine = true;
        }
    }
    
    // Debounced spectrum update to improve performance
    const updateSpectrumDebounced = debounce(function(rt) {
        // Find the closest spectrum to this retention time
        const spectrum = findClosestSpectrum(rt);
        
        // Only update if we found a spectrum and it's different from the last one
        if (spectrum && (!lastSpectrum || spectrum.rt !== lastSpectrum.rt)) {
            createSpectrumPlot(spectrum);
            lastSpectrum = spectrum;
            
            // Add RT indicator text
            const mzDiv = document.getElementById('mz-spectrum-display');
            const rtIndicator = document.createElement('div');
            rtIndicator.innerHTML = 
                `<div style="font-size:12px; color:#666; margin-top:-15px; text-align:center">
                    Cursor position: ${rt.toFixed(2)} min | 
                    Showing spectrum at RT=${spectrum.rt.toFixed(2)} min
                </div>`;
            
            // Insert as second child (after the title)
            if (mzDiv.children.length > 1) {
                mzDiv.insertBefore(rtIndicator, mzDiv.children[1]);
            } else {
                mzDiv.appendChild(rtIndicator);
            }
        }
    }, 100); // Update at most every 100ms for performance
    
    // Replace click handler with hover handler
    myPlot.on('plotly_hover', function(data) {
        if (data.points.length > 0) {
            const point = data.points[0];
            const rt = point.x;
            
            // Update the crosshair position
            updateCrosshair(rt);
            
            // Update the spectrum display
            updateSpectrumDebounced(rt);
        }
    });
    
    // Optional: Add a click handler that "locks" the current spectrum
    myPlot.on('plotly_click', function(data) {
        if (data.points.length > 0) {
            const point = data.points[0];
            const rt = point.x;
            let compoundIndex = -1;
            
            // Check for compound index in customdata
            if (point.customdata && point.customdata.length > 1) {
                compoundIndex = point.customdata[1];
            }
            
            // Find spectrum
            let spectrum = null;
            if (compoundIndex >= 0) {
                spectrum = massSpectra.find(s => s.index === compoundIndex);
            } else {
                spectrum = findClosestSpectrum(rt);
            }
            
            // Display with "locked" indicator if found
            if (spectrum) {
                createSpectrumPlot(spectrum);
                lastSpectrum = spectrum;
                
                // Add "locked" indicator
                const mzDiv = document.getElementById('mz-spectrum-display');
                const lockIndicator = document.createElement('div');
                lockIndicator.innerHTML = 
                    `<div style="font-size:12px; color:#666; margin-top:-15px; text-align:center">
                        Spectrum locked at RT=${spectrum.rt.toFixed(2)} min
                    </div>`;
                
                // Insert as second child (after the title)
                if (mzDiv.children.length > 1) {
                    mzDiv.insertBefore(lockIndicator, mzDiv.children[1]);
                } else {
                    mzDiv.appendChild(lockIndicator);
                }
            }
        }
    });
});