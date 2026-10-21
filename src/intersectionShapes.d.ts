interface InteractiveAtom {
    clickable?: boolean;
    hoverable?: boolean;
    intersectionShape?: any;
}

export function resetIntersectionShapes(atoms: InteractiveAtom[]): void;
