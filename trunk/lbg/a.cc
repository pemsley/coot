//    std::cout << "model: " << model << std::endl;
//    polyline1 = goo_canvas_polyline_new(root, TRUE, 0,
// 				       "line-width", 4.0,
// 				       "points", points,
// 				       "stroke_color", "black",
// 				       NULL);
  GooCanvasItem *rect_item, *text_item;

  rect_item = goo_canvas_rect_new (root, 100, 100, 400, 400,
				   "line-width", 10.0,
				   "radius-x", 20.0,
				   "radius-y", 10.0,
				   "stroke-color", "yellow",
				   "fill-color", "red",
				   NULL);

  text_item = goo_canvas_text_new (root, "Hello World", 300, 300, -1,
				   GTK_ANCHOR_CENTER,
				   "font", "Sans 24",
				   NULL);
  goo_canvas_item_rotate (text_item, 45, 300, 300);
  

      std::cout << "(" << atom_pos_index[i].first.x << ","
		<< atom_pos_index[i].first.y << ")"
		<< " to "
		<< "(" << atom_pos_index[j].first.x << ","
		<< atom_pos_index[j].first.y << ")"
		<< std::endl;
