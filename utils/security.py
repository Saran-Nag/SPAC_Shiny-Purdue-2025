"""
Security utilities for SPAC Shiny application.

This module provides security enhancements including HTTP security headers
and cookie security policies to mitigate various web vulnerabilities.
"""

from shiny import ui
from typing import Dict, Any


def create_security_headers() -> ui.TagList:
    """
    Create security headers to enhance application security.
    
    This function creates HTML meta tags and script elements to implement
    security policies that help protect against various attacks including
    XSS, clickjacking, and content type sniffing.
    
    Returns
    -------
    shiny.ui.TagList
        A TagList containing security-related HTML elements
        
    Notes
    -----
    These headers provide defense-in-depth security measures:
    - Content Security Policy (CSP) to prevent XSS attacks
    - X-Content-Type-Options to prevent MIME type sniffing
    - X-Frame-Options to prevent clickjacking
    - Referrer-Policy to control referrer information
    - Permissions-Policy to restrict browser features
    """
    return ui.TagList(
        # Content Security Policy to prevent XSS
        ui.tags.meta(
            **{
                "http-equiv": "Content-Security-Policy",
                "content": (
                    "default-src 'self'; "
                    "script-src 'self' 'unsafe-inline' 'unsafe-eval' "
                    "https://cdn.jsdelivr.net https://unpkg.com; "
                    "style-src 'self' 'unsafe-inline' "
                    "https://cdn.jsdelivr.net https://fonts.googleapis.com; "
                    "font-src 'self' https://fonts.gstatic.com; "
                    "img-src 'self' data: https:; "
                    "connect-src 'self'; "
                    "frame-ancestors 'none';"
                )
            }
        ),
        
        # Prevent MIME type sniffing
        ui.tags.meta(
            **{
                "http-equiv": "X-Content-Type-Options",
                "content": "nosniff"
            }
        ),
        
        # Prevent clickjacking
        ui.tags.meta(
            **{
                "http-equiv": "X-Frame-Options",
                "content": "DENY"
            }
        ),
        
        # Control referrer policy
        ui.tags.meta(
            **{
                "http-equiv": "Referrer-Policy",
                "content": "strict-origin-when-cross-origin"
            }
        ),
        
        # Restrict browser features
        ui.tags.meta(
            **{
                "http-equiv": "Permissions-Policy",
                "content": (
                    "geolocation=(), "
                    "microphone=(), "
                    "camera=(), "
                    "payment=(), "
                    "usb=(), "
                    "magnetometer=(), "
                    "gyroscope=(), "
                    "accelerometer=()"
                )
            }
        )
    )


def create_cookie_security_script() -> ui.Tag:
    """
    Create JavaScript to enhance cookie security.
    
    This function creates a script that monitors and enhances cookie security
    by attempting to modify insecure cookies to include security attributes.
    While this cannot modify HttpOnly cookies (by design), it can help with
    other security attributes.
    
    Returns
    -------
    shiny.ui.Tag
        A script tag containing JavaScript for cookie security enhancement
        
    Notes
    -----
    This is a client-side mitigation that:
    1. Monitors for insecure cookies
    2. Attempts to enhance cookie security attributes where possible
    3. Logs security warnings for monitoring purposes
    
    Limitations:
    - Cannot modify HttpOnly cookies (this is intentional browser security)
    - Cannot modify cookies from different domains (CORS policy)
    - AWS ALB cookies are typically HttpOnly and cross-domain managed
    """
    cookie_security_js = """
    (function() {
        'use strict';
        
        // Function to check and report cookie security status
        function auditCookieSecurity() {
            const cookies = document.cookie.split(';');
            const insecureCookies = [];
            
            cookies.forEach(function(cookie) {
                const cookieName = cookie.split('=')[0].trim();
                if (cookieName && !cookieName.startsWith('__Secure-')) {
                    // Log cookies that might not have security attributes
                    insecureCookies.push(cookieName);
                }
            });
            
            if (insecureCookies.length > 0) {
                console.warn('SPAC Security Audit: Found cookies without ' +
                           'secure prefix:', insecureCookies);
                console.info('Note: HttpOnly cookies (which are secure) ' +
                           'cannot be read by JavaScript');
            }
        }
        
        // Function to set secure defaults for future cookies
        function setSecureCookieDefaults() {
            // Override document.cookie setter to add security attributes
            const originalCookieDesc = 
                Object.getOwnPropertyDescriptor(Document.prototype, 'cookie') ||
                Object.getOwnPropertyDescriptor(HTMLDocument.prototype, 
                                              'cookie');
            
            if (originalCookieDesc && originalCookieDesc.set) {
                Object.defineProperty(document, 'cookie', {
                    set: function(value) {
                        // Add security attributes if not already present
                        let secureValue = value;
                        
                        // Add Secure flag for HTTPS connections
                        if (!secureValue.toLowerCase().includes('secure') && 
                            location.protocol === 'https:') {
                            secureValue += '; Secure';
                        }
                        
                        // Add SameSite attribute if not present
                        if (!secureValue.toLowerCase().includes('samesite')) {
                            secureValue += '; SameSite=Strict';
                        }
                        
                        return originalCookieDesc.set.call(this, secureValue);
                    },
                    get: originalCookieDesc.get
                });
            }
        }
        
        // Run security audit on page load
        document.addEventListener('DOMContentLoaded', function() {
            auditCookieSecurity();
            setSecureCookieDefaults();
        });
        
        // Monitor for changes (e.g., AJAX requests that might set cookies)
        setInterval(auditCookieSecurity, 30000); // Check every 30 seconds
        
    })();
    """
    
    return ui.tags.script(cookie_security_js)


def create_security_recommendations_comment() -> ui.Tag:
    """
    Create HTML comment with security recommendations for infrastructure team.
    
    Returns
    -------
    shiny.ui.Tag
        HTML comment containing security recommendations
    """
    recommendations = """
    SECURITY RECOMMENDATIONS FOR INFRASTRUCTURE TEAM:
    
    1. AWS Application Load Balancer (ALB) Cookie Security:
       - Configure ALB to set HttpOnly flag on AWSALBTG and 
         AWSALBTGCORS cookies
       - Configure ALB to set Secure flag on all cookies for HTTPS
       - This requires ALB configuration changes, not application-level
       - Reference: AWS ELB Target Group Cookie Settings Documentation
    
    2. Posit Connect (RSConnect) Cookie Security:
       - Configure rscid and rscid-legacy cookies with Secure flag
       - Ensure all session cookies have both HttpOnly and Secure flags
       - Review Posit Connect server configuration for cookie security
    
    3. Additional Security Headers (if not set at load balancer level):
       - Strict-Transport-Security: max-age=31536000; includeSubDomains
       - X-Content-Type-Options: nosniff
       - X-Frame-Options: DENY
       - Content-Security-Policy: [see application security headers]
    
    4. Cookie Security Best Practices:
       - All cookies should have HttpOnly flag when they don't need 
         JavaScript access
       - All cookies should have Secure flag when served over HTTPS
       - Use SameSite attribute to prevent CSRF attacks
       - Consider using __Secure- prefix for security-critical cookies
    
    5. Monitoring:
       - Monitor security headers using tools like securityheaders.com
       - Regular security scans with tools like OWASP ZAP or Burp Suite
       - Log and monitor cookie security violations
       - Set up alerts for insecure cookie configurations
    """
    
    return ui.HTML(f"<!-- {recommendations} -->")


def apply_security_enhancements() -> ui.TagList:
    """
    Apply comprehensive security enhancements to the application.
    
    This function combines all security measures into a single TagList
    that can be easily included in the application UI.
    
    Returns
    -------
    shiny.ui.TagList
        Complete set of security enhancements for the application
        
    Examples
    --------
    >>> from utils.security import apply_security_enhancements
    >>> 
    >>> app_ui = ui.page_fluid(
    ...     apply_security_enhancements(),
    ...     # ... rest of application UI
    ... )
    """
    return ui.TagList(
        create_security_recommendations_comment(),
        create_security_headers(),
        create_cookie_security_script()
    )


def get_security_recommendations() -> Dict[str, Any]:
    """
    Get structured security recommendations for documentation and deployment.
    
    Returns
    -------
    Dict[str, Any]
        Dictionary containing categorized security recommendations
    """
    return {
        "infrastructure": {
            "aws_alb": {
                "description": "Configure AWS Application Load Balancer cookie security",
                "actions": [
                    "Set HttpOnly flag on AWSALBTG and AWSALBTGCORS cookies",
                    "Set Secure flag on all ALB cookies for HTTPS",
                    "Configure ALB target group cookie settings",
                    "Review load balancer security group rules"
                ],
                "priority": "HIGH",
                "impact": "Prevents XSS-based session hijacking and " +
                         "man-in-the-middle cookie theft"
            },
            "posit_connect": {
                "description": "Configure Posit Connect cookie security",
                "actions": [
                    "Set Secure flag on rscid and rscid-legacy cookies",
                    "Ensure HttpOnly flag is maintained on session cookies",
                    "Review Posit Connect server cookie configuration",
                    "Update server.conf with secure cookie settings"
                ],
                "priority": "HIGH", 
                "impact": "Prevents session hijacking over insecure connections"
            },
            "security_headers": {
                "description": "Implement security headers at infrastructure level",
                "actions": [
                    "Configure HSTS header with max-age=31536000",
                    "Set X-Content-Type-Options: nosniff",
                    "Configure Content-Security-Policy",
                    "Set Referrer-Policy: strict-origin-when-cross-origin"
                ],
                "priority": "MEDIUM",
                "impact": "Defense against multiple attack vectors"
            }
        },
        "application": {
            "cookie_security": {
                "description": "Application-level cookie security enhancements",
                "actions": [
                    "Implement cookie security monitoring",
                    "Add client-side security audit logging",
                    "Override cookie defaults for security attributes",
                    "Use __Secure- prefix for critical cookies"
                ],
                "priority": "LOW",
                "impact": "Additional layer of defense and monitoring"
            },
            "content_security": {
                "description": "Content Security Policy implementation",
                "actions": [
                    "Define restrictive CSP for scripts and styles",
                    "Implement nonce-based CSP for inline scripts",
                    "Monitor CSP violations",
                    "Regularly review and update CSP directives"
                ],
                "priority": "MEDIUM",
                "impact": "Prevents XSS and code injection attacks"
            }
        },
        "monitoring": {
            "security_scanning": {
                "description": "Regular security assessments",
                "actions": [
                    "Automated security header scanning",
                    "Cookie security audit logging",
                    "Vulnerability scanning with OWASP ZAP",
                    "Penetration testing schedule",
                    "SSL/TLS configuration testing"
                ],
                "priority": "MEDIUM",
                "impact": "Early detection of security issues"
            }
        }
    }


def create_deployment_checklist() -> Dict[str, Any]:
    """
    Create a deployment security checklist for operations teams.
    
    Returns
    -------
    Dict[str, Any]
        Structured checklist for secure deployment
    """
    return {
        "pre_deployment": {
            "infrastructure": [
                "Verify HTTPS is enforced on all endpoints",
                "Configure ALB/Load Balancer cookie security settings",
                "Set up security headers at infrastructure level",
                "Configure HSTS with appropriate max-age",
                "Review and apply security group rules"
            ],
            "application": [
                "Update application with security enhancements",
                "Test cookie security in staging environment", 
                "Verify CSP doesn't break application functionality",
                "Test application with security scanning tools"
            ]
        },
        "post_deployment": {
            "verification": [
                "Run security header scan (securityheaders.com)",
                "Verify all cookies have appropriate security flags",
                "Test with OWASP ZAP or similar security scanner",
                "Check browser developer tools for CSP violations",
                "Verify HSTS is working correctly"
            ],
            "monitoring": [
                "Set up security monitoring alerts",
                "Configure log analysis for security events",
                "Schedule regular security scans",
                "Monitor for new security vulnerabilities"
            ]
        }
    }