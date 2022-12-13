//Alignment for Sprites

import { Vector2 } from "../math";

export const SpriteAlignment = {
    topLeft: new Vector2(1, -1),
    topCenter: new Vector2(0, -1),
    topRight: new Vector2(-1, -1),
    centerLeft: new Vector2(1, 0),
    center: new Vector2(0, 0),
    centerRight: new Vector2(-1, 0),
    bottomLeft: new Vector2(1, 1),
    bottomCenter: new Vector2(0, 1),
    bottomRight: new Vector2(-1, 1)
}
