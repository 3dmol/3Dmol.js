// Event Handling
export default class EventDispatcher {
  listeners = {};

  addEventListener(type, listener) {
    if (this.listeners[type] === undefined) this.listeners[type] = [];

    if (this.listeners[type].indexOf(listener) === -1) this.listeners[type].push(listener);
  }

  removeEventListener(type, listener) {
    const index = this.listeners[type].indexOf(listener);

    if (index !== -1) this.listeners[type].splice(index, 1);
  }

  dispatchEvent(event) {
    const listenerArray = this.listeners[event.type];

    if (listenerArray !== undefined) {
      event.target = this;

      for (let i = 0, l = listenerArray.length; i < l; i++) listenerArray[i].call(this, event);
    }
  }
}
