"""
Accessibility utilities for SPAC Shiny application.

This module provides common functions to ensure UI components
meet accessibility standards (WCAG, Section 508, ANDI compliance).
"""

from shiny import ui


def accessible_slider(slider_id, label, min_val, max_val, value=0, step=1,
                      aria_description=None):
    """
    Create an accessible slider with proper ARIA labeling.
    
    This function creates a slider that automatically adds accessibility
    attributes to the dynamically generated IonRangeSlider elements,
    ensuring ANDI compliance and screen reader compatibility.
    
    Args:
        slider_id (str): Unique ID for the slider input
        label (str): Visible label text for the slider
        min_val (int/float): Minimum slider value
        max_val (int/float): Maximum slider value
        value (int/float, optional): Initial slider value. Defaults to 0.
        step (int/float, optional): Step increment. Defaults to 1.
        aria_description (str, optional): Custom ARIA label.
                                        If None, uses the label text.
    
    Returns:
        ui.div: A div containing the slider with accessibility script
        
    Example:
        accessible_slider(
            "rotation_slider",
            "Rotate X-axis Labels (degrees)",
            0, 90, value=0, step=1
        )
    """
    # Generate ARIA label from the visible label if not provided
    if aria_description is None:
        aria_description = f"{label} control"
    
    return ui.div(
        ui.input_slider(
            slider_id,
            label,
            min=min_val,
            max=max_val,
            value=value,
            step=step
        ),
        _slider_accessibility_script(slider_id, aria_description)
    )


def _slider_accessibility_script(slider_id, aria_label):
    """
    Generate JavaScript to add accessibility attributes to slider.
    
    This internal function creates the JavaScript needed to find the
    dynamically generated IonRangeSlider elements and add proper
    ARIA labels for screen reader compatibility.
    
    Args:
        slider_id (str): The slider's HTML ID
        aria_label (str): The ARIA label to apply
        
    Returns:
        ui.tags.script: Script tag with accessibility JavaScript
    """
    return ui.tags.script(
        f"""
        function addAccessibilityTo_{slider_id.replace('-', '_')}() {{
            const slider = document.querySelector('#{slider_id}');
            if (!slider) return;
            
            // Find the irs-line element within the form group container
            const container = slider.closest('.form-group');
            if (!container) return;
            
            const line = container.querySelector('.irs-line[tabindex="0"]');
            if (!line) return;
            
            if (!line.hasAttribute('aria-label')) {{
                line.setAttribute('aria-label', '{aria_label}');
            }}
        }}
        
        // Multiple attempts to handle timing of IonRangeSlider initialization
        setTimeout(addAccessibilityTo_{slider_id.replace('-', '_')}, 1000);
        setTimeout(addAccessibilityTo_{slider_id.replace('-', '_')}, 2000);
        setTimeout(addAccessibilityTo_{slider_id.replace('-', '_')}, 3000);
        """
    )


def apply_slider_accessibility_global():
    """
    Apply accessibility fixes to all sliders on the page globally.
    
    This function provides a global solution that can be added once
    to fix all sliders on a page, rather than individual fixes.
    Useful for pages with many sliders.
    
    Returns:
        ui.tags.script: Global accessibility script for all sliders
    """
    return ui.tags.script(
        """
        function addGlobalSliderAccessibility() {
            // Find all IRS line elements (the actual interactive slider elements)
            const sliders = document.querySelectorAll('.irs-line[tabindex="0"]');
            
            sliders.forEach(line => {
                if (line.hasAttribute('aria-label')) {
                    return; // Already has aria-label
                }
                
                // Find the associated input to get label info
                const container = line.closest('.form-group');
                if (!container) {
                    // Try alternative approach - look for parent with slider ID
                    const parentDiv = line.closest('div[id$="_slider"], div[id*="slider"]');
                    if (parentDiv) {
                        line.setAttribute('aria-label', parentDiv.id.replace('_', ' ') + ' control');
                    } else {
                        line.setAttribute('aria-label', 'Slider control');
                    }
                    return;
                }
                
                const label = container.querySelector('label');
                const labelText = label?.textContent?.trim() || 'Slider control';
                line.setAttribute('aria-label', labelText);
            });
        }
        
        // Set up mutation observer for slider elements
        function setupSliderObserver() {
            const observer = new MutationObserver(function(mutations) {
                mutations.forEach(function(mutation) {
                    if (mutation.type === 'childList') {
                        // Check if new IRS elements were added
                        mutation.addedNodes.forEach(function(node) {
                            if (node.nodeType === 1) { // Element node
                                const newSliders = node.querySelectorAll ? 
                                    node.querySelectorAll('.irs-line[tabindex="0"]') : [];
                                if (newSliders.length > 0) {
                                    setTimeout(addGlobalSliderAccessibility, 100);
                                }
                            }
                        });
                    }
                });
            });
            
            observer.observe(document.body, {
                childList: true,
                subtree: true
            });
        }
        
        // Apply to all sliders with multiple timing attempts
        setTimeout(addGlobalSliderAccessibility, 1000);
        setTimeout(addGlobalSliderAccessibility, 2000);
        setTimeout(addGlobalSliderAccessibility, 3000);
        setTimeout(addGlobalSliderAccessibility, 5000);
        
        // Set up slider monitoring
        setTimeout(setupSliderObserver, 1000);
        """
    )


def debug_all_accessibility():
    """
    Add comprehensive debugging for all accessibility elements on the page.
    
    This function helps troubleshoot accessibility issues by logging detailed
    information about all slider and navigation elements.
    
    Returns:
        ui.tags.script: Global debugging script for accessibility elements
    """
    return ui.tags.script(
        """
        function debugAllAccessibility() {
            console.warn('=== ACCESSIBILITY DEBUG REPORT START ===');
            console.warn('Debug function is running at:', new Date());
            
            // Debug all sliders
            console.warn('--- SLIDER DEBUG ---');
            const allSliders = document.querySelectorAll('input[type="range"]');
            console.warn('Found', allSliders.length, 'range inputs');
            
            allSliders.forEach((slider, index) => {
                console.log('Slider', index + ':', slider.id);
                const container = slider.closest('.form-group');
                if (container) {
                    const irsElements = container.querySelectorAll('.irs-line');
                    console.log('  - IRS elements found:', irsElements.length);
                    irsElements.forEach((irs, irsIndex) => {
                        console.log('    IRS', irsIndex, '- tabindex:', irs.getAttribute('tabindex'), 
                                  'aria-label:', irs.getAttribute('aria-label'));
                    });
                } else {
                    console.log('  - No form-group container found');
                }
            });
            
            // Debug all IRS elements directly
            console.log('\\n--- DIRECT IRS DEBUG ---');
            const allIrsLines = document.querySelectorAll('.irs-line');
            console.log('Found', allIrsLines.length, 'total .irs-line elements');
            
            allIrsLines.forEach((line, index) => {
                console.log('IRS Line', index + ':');
                console.log('  - tabindex:', line.getAttribute('tabindex'));
                console.log('  - aria-label:', line.getAttribute('aria-label'));
                console.log('  - classes:', line.className);
                console.log('  - parent container:', line.closest('.form-group') ? 'found' : 'not found');
            });
            
            // Debug navigation tabs
            console.log('\\n--- NAVIGATION DEBUG ---');
            const navTabs = document.querySelectorAll('.nav-tabs .nav-link, .card-header .nav-link');
            console.log('Found', navTabs.length, 'navigation tabs');
            
            navTabs.forEach((tab, index) => {
                console.log('Nav Tab', index + ':');
                console.log('  - tabindex:', tab.getAttribute('tabindex'));
                console.log('  - role:', tab.getAttribute('role'));
                console.log('  - text:', tab.textContent.trim());
            });
            
            console.warn('=== END DEBUG REPORT ===');
        }
        
        // Run debug immediately and with multiple attempts
        console.warn('ACCESSIBILITY DEBUG SCRIPT LOADED!');
        alert('Debug script loaded - check console for accessibility debug info');
        debugAllAccessibility(); // Run immediately
        setTimeout(debugAllAccessibility, 500);
        setTimeout(debugAllAccessibility, 1500);
        setTimeout(debugAllAccessibility, 3500);
        setTimeout(debugAllAccessibility, 5000);
        """
    )


def accessible_navigation():
    """
    Fix navigation tab accessibility issues for navset_card_tab components.
    
    This function addresses ANDI warnings about "Focusable element is not in 
    keyboard tab order" by ensuring navigation tabs are properly accessible
    through keyboard navigation.
    
    The issue occurs because Bootstrap/Shiny sets tabindex="-1" on inactive
    nav tabs, removing them from keyboard tab order. This fix ensures all
    tabs remain keyboard accessible while maintaining proper ARIA states.
    
    Returns:
        ui.tags.script: JavaScript to fix navigation accessibility
    """
    return ui.tags.script(
        """
        function fixNavigationAccessibility() {
            // Find all navigation tabs in navset containers
            const navTabs = document.querySelectorAll(
                '.nav-tabs .nav-link, .card-header .nav-link'
            );
            
            navTabs.forEach(tab => {
                // Ensure all tabs have explicit tabindex for keyboard accessibility
                const currentTabindex = tab.getAttribute('tabindex');
                if (currentTabindex === '-1' || currentTabindex === null || currentTabindex === '') {
                    tab.setAttribute('tabindex', '0');
                }
                
                // Ensure proper ARIA attributes for screen readers
                if (!tab.hasAttribute('role')) {
                    tab.setAttribute('role', 'tab');
                }
                
                // Mark as processed to avoid re-processing
                if (!tab.hasAttribute('data-accessibility-fixed')) {
                    tab.setAttribute('data-accessibility-fixed', 'true');
                    
                    // Add keyboard event handling for accessibility
                    tab.addEventListener('keydown', function(e) {
                        // Handle Enter and Space key activation
                        if (e.key === 'Enter' || e.key === ' ') {
                            e.preventDefault();
                            tab.click();
                        }
                        
                        // Handle arrow key navigation between tabs
                        if (e.key === 'ArrowLeft' || e.key === 'ArrowRight') {
                            const tabs = Array.from(tab.closest('.nav').querySelectorAll('.nav-link'));
                            const currentIndex = tabs.indexOf(tab);
                            let nextIndex;
                            
                            if (e.key === 'ArrowRight') {
                                nextIndex = (currentIndex + 1) % tabs.length;
                            } else {
                                nextIndex = currentIndex === 0 ? tabs.length - 1 : currentIndex - 1;
                            }
                            
                            tabs[nextIndex].focus();
                            e.preventDefault();
                        }
                    });
                }
            });
        }
        
        // Set up mutation observer to watch for DOM changes and reapply fixes
        function setupNavigationObserver() {
            const observer = new MutationObserver(function(mutations) {
                mutations.forEach(function(mutation) {
                    if (mutation.type === 'attributes' && mutation.attributeName === 'tabindex') {
                        const target = mutation.target;
                        if (target.matches('.nav-link')) {
                            const currentTabindex = target.getAttribute('tabindex');
                            if (currentTabindex === '-1' || currentTabindex === null || currentTabindex === '') {
                                target.setAttribute('tabindex', '0');
                            }
                        }
                    }
                });
            });
            
            // Start observing with more comprehensive options
            observer.observe(document.body, {
                attributes: true,
                attributeFilter: ['tabindex', 'class', 'aria-selected'],
                subtree: true,
                attributeOldValue: true
            });
            
            // Also add click event listeners to all nav tabs
            const navTabs = document.querySelectorAll('.nav-tabs .nav-link, .card-header .nav-link');
            navTabs.forEach(tab => {
                tab.addEventListener('click', function() {
                    // Fix navigation
                    setTimeout(fixNavigationAccessibility, 100);
                    setTimeout(fixNavigationAccessibility, 500);
                    // Fix sliders (new content may have loaded)
                    setTimeout(addGlobalSliderAccessibility, 200);
                    setTimeout(addGlobalSliderAccessibility, 800);
                    setTimeout(addGlobalSliderAccessibility, 1500);
                });
            });
        }
        
        // Immediate fix for existing elements
        const immediateNavTabs = document.querySelectorAll('.nav-tabs .nav-link, .card-header .nav-link');
        immediateNavTabs.forEach(tab => {
            const currentTabindex = tab.getAttribute('tabindex');
            if (currentTabindex === '-1' || currentTabindex === null || currentTabindex === '') {
                tab.setAttribute('tabindex', '0');
            }
            if (!tab.hasAttribute('role')) {
                tab.setAttribute('role', 'tab');
            }
        });
        
        // Run immediately if DOM is already ready
        if (document.readyState === 'loading') {
            document.addEventListener('DOMContentLoaded', fixNavigationAccessibility);
        } else {
            fixNavigationAccessibility();
        }
        
        // Apply navigation fixes with multiple timing attempts
        setTimeout(fixNavigationAccessibility, 50);
        setTimeout(fixNavigationAccessibility, 100);
        setTimeout(fixNavigationAccessibility, 200);
        setTimeout(fixNavigationAccessibility, 500);
        setTimeout(fixNavigationAccessibility, 1000);
        setTimeout(fixNavigationAccessibility, 2000);
        
        // Set up continuous monitoring
        setTimeout(setupNavigationObserver, 500);
        """
    )
