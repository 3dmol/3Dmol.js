const editingForm = {
    isEditing: false,
}


export default new Proxy(editingForm, {
    set: (target, key, value) => {
        if (key === 'isEditing') {
            target[key] = value;
            return true;
        }
        throw new Error("editingForm.isEditing is a singleton, you can't set anything else");
    },
    get(target, key) {
        if (key === 'isEditing') {
            return target[key];
        }
        throw new Error("editingForm.isEditing is a singleton, you can't get anything else");
    }
})