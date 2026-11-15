/* Compatibility entry points for pages using the standalone viewer's old
 * classroom functions. New embeds should use StateManager's learning option or
 * LearningSession from 3Dmol.ui. Protocol, dialogs, and socket ownership now live
 * entirely in that UI module. Initialization and joining return promises.
 */
var initSessions;
var joinSession;
(function() {
    var app = null, pending = null, session = null, generation = 0;
    initSessions = function() {
        var current = window.viewerApp;
        if (!current || typeof current.enableSessions !== 'function') {
            return Promise.reject(new Error('Initialize the standalone viewer before enabling classroom sessions.'));
        }
        if (current !== app) {
            app = current; pending = null; session = null; generation++;
        }
        if (!pending) {
            var loading = Promise.resolve().then(function() {
                return current.enableSessions();
            }).then(function(value) {
                if (app === current && pending === loading) session = value;
                return value;
            }).catch(function(error) {
                if (app === current && pending === loading) pending = null;
                throw error;
            });
            pending = loading;
        }
        return pending;
    };
    joinSession = function(name) {
        var initialized = initSessions(), token = generation;
        return initialized.then(function(value) {
            if (token === generation) return value.join(name);
        });
    };
    window.leaveViewerSession = function() {
        generation++;
        if (session) return session.leave();
        // Closing while support loads also cancels a queued join.
        if (pending) return pending.then(function(value) { return value.leave(); });
    };
})();
