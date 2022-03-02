//Event Handling
/** @this {EventDispatcher} */
export class EventDispatcher {
    listeners = {};
    addEventListener(type, listener) {
        if (listeners[type] === undefined)
            listeners[type] = [];

        if (listeners[type].indexOf(listener) === -1)
            listeners[type].push(listener);
    };

    removeEventListener(type, listener) {

        var index = listeners[type].indexOf(listener);

        if (index !== -1)
            listeners[type].splice(index, 1);

    };

    dispatchEvent(event) {

        var listenerArray = listeners[event.type];

        if (listenerArray !== undefined) {
            event.target = this;

            for (var i = 0, l = listenerArray.length; i < l; i++)
                listenerArray[i].call(this, event);

        }

    }
};